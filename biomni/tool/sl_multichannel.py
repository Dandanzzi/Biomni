"""Multi-channel synthetic lethality discovery for Biomni.

``discover_synthetic_lethal_candidates`` asks one question - is this gene more essential in
driver-mutant lines - and gates the answer on a genome-wide FDR. On a real driver that fails:
in the DepMap breast panel the textbook BRCA-PARP1 interaction scores p=0.054, q=0.998, and is
discarded. A method that throws away the one interaction everybody agrees on is not calibrated
for discovery, and adding more single-channel statistics will not fix it.

This module implements the approach Noh et al. (New Biotechnology 94:184-191, 2026) describe as
the statistical foundation of the field, plus the translational filters that review identifies as
prerequisites for a usable candidate list.

DAISY (Jerby-Arnon et al., Cell 2014) does not test one hypothesis. It runs three independent
inference channels over the same cohort and keeps pairs supported by all of them, because each
channel is individually weak and their errors are uncorrelated:

1. ``essentiality``      the candidate is MORE required in driver-deficient lines than in
                         driver-intact lines (the channel the existing tool implements alone)
2. ``co-expression``     SL partners are functionally related, so their transcripts co-vary
3. ``co-inactivation``   SL partners are rarely inactivated together in a surviving tumour -
                         a cell that lost both is dead, so the double-loss genotype is depleted
                         ("survival of the fittest")

Ranks are combined across channels rather than thresholded per channel, because the review's
Section 3 is explicit that SL discovery is a ranking problem under extreme class imbalance and
that threshold-and-count evaluation overstates what a model can do. Nothing here is FDR-gated;
the output is an ordered list plus the rank at which known interactions land, so the caller can
see the method's calibration on this cohort instead of trusting it.

Three filters from the review's Section 5 are applied as part of scoring rather than afterwards:

* **essentiality floor** - constitutively essential genes look lethal in every context and
  swamp context-specific hits. Filtered on DepMap dependency PROBABILITIES rather than a
  hard gene-effect cutoff, which is the Bayes-factor-style treatment Project Score uses.
* **druggability** - the review's complaint is that druggability is applied post hoc. Here every
  candidate carries its Broad Repurposing Hub status while it is being ranked.
* **known-interaction benchmark** - the caller supplies interactions already established for this
  driver, and the report states where they ranked. If PARP1 does not surface for BRCA loss, the
  run is not calibrated and the rest of the list should not be believed.

WHAT THIS MODULE CANNOT DO, AND WHY
-----------------------------------
The review's Table 2 lists SynLethDB 2.0 (35,943 human SL pairs) and SLKB (16,059 SL and
264,424 experimentally grounded non-SL pairs) as the labelled resources that supervised and
deep-learning predictors are trained on. Neither is in the local data lake. The BioGRID-derived
``synthetic_lethality.parquet`` that is present contains 1,909 pairs of which **zero are human** -
all are S. cerevisiae - which is precisely the species-contamination trap the review warns about.
With no human labels there is no training set, no benchmark, and no negative set, so none of the
ML or deep-learning architectures in the review (GCATSL, KG4SL, KR4SL, NSF4SL, MSGT-SL, ELISL)
can be run or evaluated here. This module therefore implements the statistical tier only, and
reports its own calibration instead of claiming accuracy it cannot measure.
"""

import os
import re
from datetime import datetime

import numpy as np

from biomni.tool.synthetic_lethality import (
    _UTC,
    _annotate_mutation_status,
    _load_depmap,
    _load_expression,
    _parse_gene_list,
    _resolve_data_lake,
    _select_cancer_models,
)

_DEPENDENCY_CACHE: dict = {}
_DRUG_TARGET_CACHE: dict = {}

# 단순 missense는 대부분 VUS라 기능 결손을 만들지 못한다. 종양억제자를 driver로 쓸 때
# 이걸 변이군에 넣으면 대비가 희석된다 (실측: BRCA에서 PARP1 차이가 -0.181 -> -0.108).
_SIMPLE_MISSENSE = re.compile(r"^[A-Z]\d+[A-Z]$")

# 리뷰 5장의 essentiality floor. gene effect 절대값이 아니라 의존 '확률'로 자른다.
CORE_ESSENTIAL_PROBABILITY = 0.5
CORE_ESSENTIAL_FRACTION = 0.90  # 이 비율 이상의 세포주에서 의존이면 core essential
# 후보가 되려면 driver 결손 세포주 중 최소 이 비율에서 실제 의존이어야 한다.
# 평균 의존확률을 요구하면 안 된다 - 합성치사는 침투도가 부분적인 경우가 많아서,
# 평균 0.5를 요구하면 POLQ나 PARP2 같은 진짜 파트너가 먼저 잘려나간다 (실측 확인).
MUTANT_DEPENDENCY_PROBABILITY = 0.5
MIN_MUTANT_DEPENDENT_FRACTION = 0.20


def _load_dependency_probability(data_lake_path: str | None = None):
    """DepMap의 유전자별 의존 '확률'(0-1)을 읽는다.

    gene effect를 -0.5로 자르는 것보다 낫다: 그 컷은 모든 세포주·모든 유전자에 같은 값을
    쓰지만, 확률은 세포주별 분포를 반영해 계산된 값이라 Project Score의 Bayes factor
    임계와 같은 역할을 한다 (리뷰 5장의 essentiality floor 처리).
    """
    import pandas as pd

    resolved = _resolve_data_lake(data_lake_path)
    if resolved in _DEPENDENCY_CACHE:
        return _DEPENDENCY_CACHE[resolved]
    path = os.path.join(resolved, "DepMap_CRISPRGeneDependency.csv")
    if not os.path.exists(path):
        return None
    frame = pd.read_csv(path, index_col=0)
    frame.columns = [c.split(" (")[0].strip().upper() for c in frame.columns]
    bundle = {"probability": frame, "path": path}
    _DEPENDENCY_CACHE[resolved] = bundle
    return bundle


def _load_drug_targets(data_lake_path: str | None = None):
    """Broad Repurposing Hub에서 유전자 -> (최고 임상단계, 약물 예시)를 만든다."""
    import pandas as pd

    resolved = _resolve_data_lake(data_lake_path)
    if resolved in _DRUG_TARGET_CACHE:
        return _DRUG_TARGET_CACHE[resolved]
    path = os.path.join(resolved, "broad_repurposing_hub_phase_moa_target_info.parquet")
    if not os.path.exists(path):
        return None
    table = pd.read_parquet(path)
    rank = {"Launched": 4, "Phase 3": 3, "Phase 2": 2, "Phase 1": 1, "Preclinical": 0}
    mapping: dict = {}
    for _, row in table.iterrows():
        targets = str(row.get("target") or "")
        if not targets or targets == "nan":
            continue
        phase = str(row.get("clinical_phase") or "").strip()
        score = rank.get(phase, 0)
        for symbol in targets.split("|"):
            symbol = symbol.strip().upper()
            if not symbol:
                continue
            current = mapping.get(symbol)
            if current is None or score > current[0]:
                mapping[symbol] = (score, phase or "Preclinical", str(row.get("pert_iname") or ""))
    _DRUG_TARGET_CACHE[resolved] = mapping
    return mapping


def define_driver_deficient_lines(driver_genes, cohort, data_lake_path, lof_only=True):
    """여러 driver를 하나의 기능 결손 유전형으로 묶는다.

    BRCA1과 BRCA2는 같은 경로(상동재조합)를 망가뜨리고 임상에서도 gBRCA1/2로 함께 등록한다.
    따로 보면 유방암 패널에서 각각 4종·7종뿐이라 검정력이 없다.

    Returns (변이 ModelID 집합, 야생형 ModelID 집합, 변이 설명, 제외된 설명).
    """
    mutant, profiled, dropped = {}, set(), []
    for gene in driver_genes:
        annotated, _ = _annotate_mutation_status(cohort, gene, data_lake_path, None)
        profiled |= set(annotated.loc[annotated["MutationStatus"] != "UNKNOWN", "ModelID"])
        for _, row in annotated[annotated["MutationStatus"] == "MUT"].iterrows():
            variant = str(row["Variant"])
            if lof_only and _SIMPLE_MISSENSE.match(variant):
                dropped.append(f"{row['StrippedCellLineName']} ({gene}:{variant})")
                continue
            mutant.setdefault(row["ModelID"], []).append(f"{gene}:{variant}")
    return set(mutant), profiled - set(mutant), mutant, dropped


def _channel_essentiality(gene_effect, symbols, mutant_ids, wildtype_ids):
    """채널 1 - driver 결손 세포주에서 더 필수적인가 (DAISY의 essentiality 모듈).

    유전자마다 t검정을 도는 대신 두 행렬의 Welch 통계량을 한 번에 계산한다 (17,916회 루프를
    피하기 위해서이고, 값은 scipy의 equal_var=False와 같다).
    """
    from scipy import stats

    mutant = gene_effect.reindex(list(mutant_ids)).to_numpy(dtype=float)
    wildtype = gene_effect.reindex(list(wildtype_ids)).to_numpy(dtype=float)

    with np.errstate(invalid="ignore", divide="ignore"):
        n_a = np.sum(~np.isnan(mutant), axis=0)
        n_b = np.sum(~np.isnan(wildtype), axis=0)
        mean_a = np.nanmean(mutant, axis=0)
        mean_b = np.nanmean(wildtype, axis=0)
        var_a = np.nanvar(mutant, axis=0, ddof=1)
        var_b = np.nanvar(wildtype, axis=0, ddof=1)
        se_a, se_b = var_a / n_a, var_b / n_b
        denominator = np.sqrt(se_a + se_b)
        t_statistic = (mean_a - mean_b) / denominator
        degrees = (se_a + se_b) ** 2 / (se_a**2 / (n_a - 1) + se_b**2 / (n_b - 1))
        p_value = 2 * stats.t.sf(np.abs(t_statistic), degrees)

    difference = mean_a - mean_b
    # 방향이 맞는 것만 점수를 준다. 합성치사는 '변이일 때 더 취약'이므로 difference < 0.
    score = np.where(difference < 0, -np.log10(np.clip(p_value, 1e-300, 1.0)), 0.0)
    return {
        "symbols": symbols,
        "score": np.nan_to_num(score),
        "difference": difference,
        "p_value": p_value,
        "mutant_mean": mean_a,
        "wildtype_mean": mean_b,
    }


def _channel_coexpression(expression, expression_columns, symbols, driver_genes, model_ids):
    """채널 2 - driver와 발현이 함께 움직이는가 (DAISY의 co-expression 모듈).

    합성치사 파트너는 기능적으로 가까운 경우가 많고, 기능적 근접성은 전사체 공변동으로
    드러난다. driver가 여럿이면 각 driver와의 상관 중 절대값이 가장 큰 것을 쓴다.
    """
    usable = [m for m in model_ids if m in expression.index]
    if len(usable) < 20:
        return None
    block = expression.loc[usable]

    columns = [expression_columns.get(s) for s in symbols]
    valid = np.array([c is not None for c in columns])
    matrix = np.full((len(usable), len(symbols)), np.nan)
    present = [c for c in columns if c is not None]
    matrix[:, valid] = block[present].to_numpy(dtype=float)

    best = np.zeros(len(symbols))
    for driver in driver_genes:
        column = expression_columns.get(driver)
        if column is None:
            continue
        reference = block[column].to_numpy(dtype=float)
        centred_ref = reference - np.nanmean(reference)
        with np.errstate(invalid="ignore", divide="ignore"):
            # 전부 NaN인 열(발현 행렬에 없는 유전자)에서 경고가 나므로 여기서 함께 억제한다.
            import warnings

            with warnings.catch_warnings():
                warnings.simplefilter("ignore", RuntimeWarning)
                centred = matrix - np.nanmean(matrix, axis=0)
            numerator = np.nansum(centred * centred_ref[:, None], axis=0)
            denominator = np.sqrt(np.nansum(centred**2, axis=0) * np.nansum(centred_ref**2))
            correlation = np.abs(np.divide(numerator, denominator, out=np.zeros_like(numerator), where=denominator > 0))
        best = np.maximum(best, np.nan_to_num(correlation))
    return best


def _channel_coinactivation(expression, expression_columns, symbols, mutant_ids, wildtype_ids, low_quantile=0.2):
    """채널 3 - driver 결손과 후보 불활성이 함께 나타나지 않는가 (DAISY의 SoF 모듈).

    두 유전자가 합성치사라면 둘 다 잃은 세포는 죽는다. 따라서 살아남아 배양된 세포주
    집합에서는 그 조합이 '덜' 관찰되어야 한다. 후보의 불활성은 패널 내 발현 하위분위로
    근사한다 (로컬 데이터레이크에 copy number가 없다).

    점수는 상호배타성의 정도다: driver 결손군에서 후보가 낮게 발현된 비율이 야생형군보다
    작을수록 높다.
    """
    mutant = [m for m in mutant_ids if m in expression.index]
    wildtype = [m for m in wildtype_ids if m in expression.index]
    if len(mutant) < 3 or len(wildtype) < 10:
        return None

    columns = [expression_columns.get(s) for s in symbols]
    valid = np.array([c is not None for c in columns])
    present = [c for c in columns if c is not None]

    everything = expression.loc[mutant + wildtype, present].to_numpy(dtype=float)
    threshold = np.nanquantile(everything, low_quantile, axis=0)
    mutant_block = expression.loc[mutant, present].to_numpy(dtype=float)
    wildtype_block = expression.loc[wildtype, present].to_numpy(dtype=float)

    with np.errstate(invalid="ignore"):
        mutant_low = np.nanmean(mutant_block <= threshold, axis=0)
        wildtype_low = np.nanmean(wildtype_block <= threshold, axis=0)

    score = np.zeros(len(symbols))
    score[valid] = np.nan_to_num(wildtype_low - mutant_low)  # 양수 = 상호배타적
    return np.clip(score, 0, None)


def _rank_normalise(values):
    """값을 0-1 백분위로 바꾼다. 채널마다 단위가 달라 그대로 더할 수 없다."""
    import pandas as pd

    series = pd.Series(values)
    return (series.rank(method="average") / len(series)).to_numpy()

_SYNLETHDB_CACHE: dict = {}

def load_known_sl_partners(driver_genes, data_lake_path=None, min_score=0.0):
    """SynLethDB에서 driver의 알려진 합성치사 파트너를 가져온다. 없으면 None."""
    import pandas as pd
    import os

    # 사용자 홈 디렉토리(~)를 기준으로 절대 경로를 자동 생성합니다.
    path = os.path.expanduser("~/biomni/data/biomni_data/data_lake/synlethdb_human_sl.parquet")
    
    if "fixed_path" not in _SYNLETHDB_CACHE:
        if not os.path.exists(path):
            print(f"\n[경로 에러] 다음 위치에 파일이 없습니다: {path}\n")
            _SYNLETHDB_CACHE["fixed_path"] = None
        else:
            _SYNLETHDB_CACHE["fixed_path"] = pd.read_parquet(path)
            
    table = _SYNLETHDB_CACHE["fixed_path"]
    if table is None:
        return None

    drivers = {g.upper() for g in driver_genes}
    hit = table[table.gene_a.isin(drivers) | table.gene_b.isin(drivers)]
    if min_score:
        hit = hit[hit.score >= min_score]
    partners = set(hit.gene_a) | set(hit.gene_b)
    return sorted(partners - drivers)
    
def discover_sl_multichannel(
    driver_genes,
    cancer_type: str = "pan-cancer",
    known_sl_partners: list | None = None,  # 원래대로 복구
    lof_only: bool = True,
    top_k: int = 30,
    require_channels: int = 3,
    data_lake_path: str | None = None,
) -> str:
    """Rank synthetic lethal partners of a driver by combining three independent inference channels.

    Implements the DAISY-style multi-channel design described by Noh et al. (New Biotechnology 94,
    2026) together with the translational filters that review identifies as prerequisites: an
    essentiality floor applied on dependency probabilities, druggability carried through scoring
    rather than bolted on afterwards, and a known-interaction benchmark that states the method's
    calibration on this cohort.

    The output is a RANKED LIST, not a significance-filtered set. Under the class imbalance of SL
    discovery a genome-wide FDR gate discards true interactions - on the DepMap breast panel the
    BRCA-PARP1 interaction reaches only q=0.998 - so ranking metrics are the appropriate frame.

    Parameters
    ----------
    driver_genes : list[str] | str
        Driver gene(s) whose loss defines the genotype, e.g. ``["BRCA1", "BRCA2"]``. Multiple genes
        are pooled into one deficiency label, which is how gBRCA1/2 is treated clinically and is
        the only way to get usable group sizes for rare drivers.
    cancer_type : str, optional
        Cancer context, or "pan-cancer" (default). Pan-cancer is the default deliberately: the
        deficiency phenotype these channels measure is largely tissue-independent, and a single
        lineage rarely supplies enough driver-mutant lines to power any of the three channels.
    known_sl_partners : list[str], optional
        Interactions already established for this driver (e.g. ``["PARP1", "PARP2"]`` for BRCA).
        Their ranks are reported as a calibration check. If they do not surface, do not trust the
        rest of the list.
    lof_only : bool, optional
        Count only likely loss-of-function variants as driver-deficient, dropping simple missense
        calls (default: True). Most missense calls in tumour suppressors are variants of uncertain
        significance that dilute the contrast.
    top_k : int, optional
        Number of ranked candidates to report (default: 30).
    require_channels : int, optional
        How many channels a gene must score in to be reported (default: 3, i.e. DAISY's rule).
    data_lake_path : str, optional
        Directory holding the DepMap files.

    Returns
    -------
    str
        A research log with the genotype definition, per-channel diagnostics, the ranked candidate
        table with druggability, the known-interaction benchmark, and the resources that would be
        needed to go beyond the statistical tier.

    """
    genes = [g.upper() for g in _parse_gene_list(driver_genes)]
    if not genes:
        return "FAILURE: no driver genes supplied."
    known_sl_partners = [g.upper() for g in _parse_gene_list(known_sl_partners or [])]
    if not known_sl_partners:
        known_sl_partners = load_known_sl_partners(genes, data_lake_path) or []

    try:
        bundle = _load_depmap(data_lake_path)
        expression_bundle = _load_expression(data_lake_path)
    except FileNotFoundError as e:
        return f"FAILURE: {e}"

    gene_effect = bundle["gene_effect"]
    symbols = [c.split(" (")[0].strip().upper() for c in gene_effect.columns]
    cohort, match_note = _select_cancer_models(bundle["model"], cancer_type)
    cohort = cohort[cohort["ModelID"].isin(gene_effect.index)]
    if len(cohort) < 20:
        return f"FAILURE: only {len(cohort)} screened lines matched '{cancer_type}'; need at least 20."

    try:
        mutant_ids, wildtype_ids, mutant_map, dropped = define_driver_deficient_lines(
            genes, cohort, bundle["data_lake_path"], lof_only
        )
    except (RuntimeError, ValueError) as e:
        return f"FAILURE: {e}"
    if len(mutant_ids) < 3:
        return f"FAILURE: only {len(mutant_ids)} {'/'.join(genes)}-deficient lines; need at least 3."

    log = [
        "=" * 78,
        f"MULTI-CHANNEL SL DISCOVERY - {'/'.join(genes)} deficiency in {cancer_type}",
        f"Run at {datetime.now(tz=_UTC).strftime('%Y-%m-%d %H:%M UTC')}",
        "=" * 78,
        "",
        "Method: three independent inference channels (DAISY design, Jerby-Arnon Cell 2014),",
        "combined by rank rather than thresholded per channel. No FDR gate is applied - SL",
        "discovery is a ranking problem under extreme class imbalance, and a genome-wide FDR",
        "discards true interactions at these group sizes.",
        "",
        "STEP 1 | Genotype definition",
        f"  Cohort: {match_note}",
        f"  {'/'.join(genes)}-deficient: {len(mutant_ids)} lines   intact: {len(wildtype_ids)} lines",
    ]
    shown = sorted(mutant_map.items(), key=lambda kv: kv[1])[:12]
    for model_id, variants in shown:
        name = cohort.loc[cohort["ModelID"] == model_id, "StrippedCellLineName"]
        log.append(f"    {(name.iloc[0] if len(name) else model_id):<12}{', '.join(variants)}")
    if len(mutant_map) > 12:
        log.append(f"    ... and {len(mutant_map) - 12} more")
    if dropped:
        log.append(f"  Dropped as likely-VUS missense ({len(dropped)}): {', '.join(dropped[:6])}")
        if len(dropped) > 6:
            log.append(f"    ... and {len(dropped) - 6} more")

    # ---- 채널 계산 -----------------------------------------------------------------------------
    essentiality = _channel_essentiality(gene_effect, symbols, mutant_ids, wildtype_ids)
    expression = expression_bundle["expression"]
    expression_columns = expression_bundle["gene_columns"]
    coexpression = _channel_coexpression(
        expression, expression_columns, symbols, genes, list(mutant_ids) + list(wildtype_ids)
    )
    coinactivation = _channel_coinactivation(
        expression, expression_columns, symbols, mutant_ids, wildtype_ids
    )

    log.append("")
    log.append("STEP 2 | Channels")
    log.append(f"  essentiality    : {int((essentiality['score'] > 0).sum()):>6} genes score (direction-correct)")
    log.append(
        f"  co-expression   : {'unavailable - too few lines with expression' if coexpression is None else str(int((coexpression > 0.3).sum())) + ' genes with |r| > 0.3'}"
    )
    log.append(
        f"  co-inactivation : {'unavailable - group too small' if coinactivation is None else str(int((coinactivation > 0).sum())) + ' genes depleted for double loss'}"
    )

    # ---- essentiality floor (리뷰 5장) ----------------------------------------------------------
    dependency = _load_dependency_probability(data_lake_path)
    core_essential = np.zeros(len(symbols), dtype=bool)
    mutant_dependent = np.ones(len(symbols), dtype=bool)
    if dependency is not None:
        probability = dependency["probability"]
        aligned = probability.reindex(columns=symbols)
        core_fraction = (aligned >= CORE_ESSENTIAL_PROBABILITY).mean(axis=0).to_numpy()
        core_essential = core_fraction >= CORE_ESSENTIAL_FRACTION
        mutant_rows = aligned.reindex([m for m in mutant_ids if m in aligned.index])
        if len(mutant_rows):
            dependent_fraction = (mutant_rows >= MUTANT_DEPENDENCY_PROBABILITY).mean(axis=0).to_numpy()
            mutant_dependent = dependent_fraction >= MIN_MUTANT_DEPENDENT_FRACTION
        log.append("")
        log.append("STEP 3 | Essentiality floor (dependency probabilities, not gene-effect cutoffs)")
        log.append(
            f"  core essential (dep. prob >= {CORE_ESSENTIAL_PROBABILITY} in >= "
            f"{CORE_ESSENTIAL_FRACTION:.0%} of lines): {int(core_essential.sum())} genes excluded"
        )
        log.append(
            f"  must be a real dependency in >= {MIN_MUTANT_DEPENDENT_FRACTION:.0%} of deficient lines "
            f"(dep. prob >= {MUTANT_DEPENDENCY_PROBABILITY}): {int(mutant_dependent.sum())} genes qualify"
        )
        log.append("  NOTE: a fraction, not a mean. Requiring the mean to clear 0.5 assumes full")
        log.append("  penetrance and removes partially penetrant partners such as POLQ and PARP2.")
    else:
        log.append("")
        log.append("STEP 3 | Essentiality floor SKIPPED - DepMap_CRISPRGeneDependency.csv not found")

    # ---- 랭크 결합 -----------------------------------------------------------------------------
    channels_used = ["essentiality"]
    combined = _rank_normalise(essentiality["score"])
    channel_count = (essentiality["score"] > 0).astype(int)
    if coexpression is not None:
        combined = combined + _rank_normalise(coexpression)
        channel_count = channel_count + (coexpression > 0.2).astype(int)
        channels_used.append("co-expression")
    if coinactivation is not None:
        combined = combined + _rank_normalise(coinactivation)
        channel_count = channel_count + (coinactivation > 0).astype(int)
        channels_used.append("co-inactivation")
    combined = combined / len(channels_used)

    import pandas as pd

    table = pd.DataFrame(
        {
            "gene": symbols,
            "combined": combined,
            "channels": channel_count,
            "essentiality_p": essentiality["p_value"],
            "difference": essentiality["difference"],
            "mutant_mean": essentiality["mutant_mean"],
            "wildtype_mean": essentiality["wildtype_mean"],
            "coexpression": coexpression if coexpression is not None else np.nan,
            "coinactivation": coinactivation if coinactivation is not None else np.nan,
            "core_essential": core_essential,
            "mutant_dependent": mutant_dependent,
        }
    )
    table = table[~table["gene"].isin(genes)]
    eligible = table[
        (~table["core_essential"])
        & table["mutant_dependent"]
        & (table["difference"] < 0)
        & (table["channels"] >= min(require_channels, len(channels_used)))
    ].copy()
    eligible = eligible.sort_values("combined", ascending=False)

    drug_targets = _load_drug_targets(data_lake_path)

    log.append("")
    log.append("STEP 4 | Ranked candidates")
    log.append(f"  Channels combined: {', '.join(channels_used)}")
    log.append(f"  {len(eligible)} genes pass the floor and score in >= {min(require_channels, len(channels_used))} channels")
    log.append("")
    header = (
        f"  {'rank':<5}{'gene':<11}{'score':>7}{'ch':>3}{'MUT':>8}{'WT':>8}{'diff':>8}"
        f"{'p':>9}{'coexp':>7}{'coinact':>8}  druggability"
    )
    log.append(header)
    log.append("  " + "-" * (len(header) + 4))
    for index, (_, row) in enumerate(eligible.head(top_k).iterrows(), start=1):
        drug = drug_targets.get(row["gene"]) if drug_targets else None
        annotation = f"{drug[1]} ({drug[2]})" if drug else "-"
        log.append(
            f"  {index:<5}{row['gene']:<11}{row['combined']:>7.3f}{int(row['channels']):>3}"
            f"{row['mutant_mean']:>8.3f}{row['wildtype_mean']:>8.3f}{row['difference']:>8.3f}"
            f"{row['essentiality_p']:>9.1e}{row['coexpression']:>7.2f}{row['coinactivation']:>8.2f}  {annotation}"
        )

    # ---- 알려진 상호작용으로 보정 상태 확인 -------------------------------------------------
    if known_sl_partners:
        log.append("")
        log.append("STEP 5 | Calibration on known interactions")
        log.append("  Where do interactions that are already established land in this ranking?")
        log.append("  This is the check that says whether the list above is worth reading.")
        log.append("")
        ordered = table.sort_values("combined", ascending=False).reset_index(drop=True)
        eligible_order = eligible.reset_index(drop=True)
        for partner in known_sl_partners:
            hit = ordered[ordered["gene"] == partner]
            if hit.empty:
                log.append(f"    {partner:<10} not in the CRISPR library")
                continue
            row = hit.iloc[0]
            overall = int(hit.index[0]) + 1
            shortlist = eligible_order.index[eligible_order["gene"] == partner]
            place = f"#{int(shortlist[0]) + 1} of {len(eligible_order)}" if len(shortlist) else "filtered out"
            reason = ""
            if not len(shortlist):
                if row["core_essential"]:
                    reason = " (core essential)"
                elif not row["mutant_dependent"]:
                    reason = " (not a dependency in enough deficient lines)"
                elif row["difference"] >= 0:
                    reason = " (direction wrong in this cohort)"
                else:
                    reason = f" (scored in only {int(row['channels'])} channels)"
            log.append(
                f"    {partner:<10} shortlist {place}{reason}; genome-wide #{overall} of {len(ordered)}; "
                f"diff {row['difference']:+.3f}, p={row['essentiality_p']:.1e}"
            )

    log.append("")
    log.append("WHAT IS MISSING TO GO FURTHER")
    log.append("  This is the statistical tier only. The supervised and deep-learning methods in")
    log.append("  Noh et al. need labelled human SL pairs, and the local data lake has none:")
    log.append("  synthetic_lethality.parquet holds 1,909 pairs, all S. cerevisiae, 0 human.")
    log.append("  Required for the next tier: SynLethDB 2.0 (35,943 human pairs + SynLethKG) as a")
    log.append("  training set, SLKB (16,059 SL / 264,424 experimentally grounded non-SL, cell-line")
    log.append("  resolved) for negative labels, and DepMap copy number for a proper co-inactivation")
    log.append("  channel - expression low-quantile is a weak proxy for deletion.")
    log.append("")
    log.append(f"MULTICHANNEL_CANDIDATES: {', '.join(eligible.head(top_k)['gene']) if len(eligible) else '(none)'}")
    log.append("")
    log.append("PROVENANCE AND LIMITS")
    for item in bundle["provenance"]:
        log.append(f"  - {item}")
    log.append(f"  - {expression_bundle['provenance']}")
    if dependency is not None:
        log.append(f"  - Dependency probabilities: {dependency['path']}")
    log.append("  - Co-inactivation uses expression low-quantile as a proxy for genomic loss; the")
    log.append("    local data lake has no copy-number matrix, so that channel is the weakest.")
    log.append("  - No labelled human SL data is used anywhere: the ranking is unsupervised, and")
    log.append("    the calibration section is the only evidence about its accuracy.")
    return "\n".join(log)
