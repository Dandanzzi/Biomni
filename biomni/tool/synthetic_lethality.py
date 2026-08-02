"""Synthetic lethality (SL) discovery and evidence-validation tools for Biomni.

This module implements a three-stage tool chain for context-specific synthetic
lethality research (e.g. "KRAS-mutant pancreatic cancer"):

1. ``discover_synthetic_lethal_candidates`` - statistical discovery of candidate
   SL partners from DepMap CRISPR gene-effect data stratified by mutation status.
2. ``validate_sl_candidates_with_pubmed`` - literature validation of candidates via
   the NCBI Entrez E-utilities API (esearch + efetch abstract parsing).
3. ``analyze_ppi_network_for_sl`` - protein-protein interaction validation via the
   STRING-DB REST API.
4. ``check_dependency_confounders`` - tests whether a candidate's dependency is really
   explained by the driver mutation or by a confounder (paralog loss, lineage, expression
   biomarker), before any wet-lab commitment.
5. ``generate_sl_evidence_dossier`` - orchestrates 1-4 into a single evidence dossier
   with Go / Hold / No-go recommendations.

Design notes
------------
* All numeric judgements are made by explicit code (statistics, thresholds, QC),
  never by free-text summarisation, so the LLM layer only plans and narrates.
* Every tool reports provenance (data file, release, sample counts, API endpoint)
  and QC warnings (pan-essentiality, low sample size, lineage confounding) so that
  weak candidates can be actively falsified rather than silently promoted.
"""

import json
import os
import re
import time
import urllib.parse
import xml.etree.ElementTree as ET
from datetime import datetime, timezone

import requests

# Aliased so the module also imports on Python 3.10; datetime.UTC only exists from 3.11.
_UTC = timezone.utc  # noqa: UP017

DEFAULT_DATA_LAKE = "./data/biomni_data/data_lake"
CBIOPORTAL_API = "https://www.cbioportal.org/api"
CCLE_STUDY_ID = "ccle_broad_2019"
STRING_API = "https://string-db.org/api"
ENTREZ_API = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"

# Cell lines whose gene effect is below this value are treated as "depleted".
DEPLETION_THRESHOLD = -0.5
# A gene depleted in more than this fraction of ALL screened lines is pan-essential.
PAN_ESSENTIAL_FRACTION = 0.80

# Module level caches so that the ~400 MB DepMap matrices are read only once per session.
_DEPMAP_CACHE: dict = {}
_MUTATION_CACHE: dict = {}
_EXPRESSION_CACHE: dict = {}


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------
def _normalize_cell_line_name(name: str) -> str:
    """Normalise a cell line name for cross-resource matching (MIA PaCa-2 -> MIAPACA2)."""
    return re.sub(r"[^A-Z0-9]", "", str(name).upper())


def _http_get(url: str, params: dict | None = None, timeout: int = 60, retries: int = 3, as_json: bool = True):
    """GET with exponential backoff. Returns parsed JSON, raw text, or None on failure."""
    last_error = None
    for attempt in range(retries):
        try:
            response = requests.get(url, params=params, timeout=timeout)
            response.raise_for_status()
            return response.json() if as_json else response.text
        except Exception as e:  # network layer: any failure is retried
            last_error = e
            time.sleep(1.5 * (attempt + 1))
    print(f"[synthetic_lethality] request failed after {retries} attempts ({url}): {last_error}")
    return None


def _resolve_data_lake(data_lake_path: str | None) -> str:
    """Return a usable data lake directory containing the DepMap files."""
    candidates = [
        data_lake_path,
        DEFAULT_DATA_LAKE,
        os.path.join(os.path.expanduser("~"), "biomni_data", "data_lake"),
        "./biomni_data/data_lake",
    ]
    for candidate in candidates:
        if candidate and os.path.exists(os.path.join(candidate, "DepMap_CRISPRGeneEffect.csv")):
            return candidate
    raise FileNotFoundError(
        "DepMap_CRISPRGeneEffect.csv / DepMap_Model.csv were not found. Provide `data_lake_path` "
        f"explicitly (searched: {[c for c in candidates if c]})."
    )


def _file_provenance(path: str) -> str:
    stat = os.stat(path)
    modified = datetime.fromtimestamp(stat.st_mtime, tz=_UTC).strftime("%Y-%m-%d")
    return f"{os.path.basename(path)} ({stat.st_size / 1e6:.0f} MB, local snapshot dated {modified})"


def _load_depmap(data_lake_path: str | None = None) -> dict:
    """Load and cache the DepMap CRISPR gene-effect matrix and the model annotation table.

    Returns
    -------
    dict with keys ``gene_effect`` (DataFrame: ModelID x gene), ``model`` (DataFrame),
    ``gene_columns`` (mapping HUGO symbol -> matrix column) and ``provenance`` (list of str).
    """
    import pandas as pd

    resolved = _resolve_data_lake(data_lake_path)
    if resolved in _DEPMAP_CACHE:
        return _DEPMAP_CACHE[resolved]

    effect_path = os.path.join(resolved, "DepMap_CRISPRGeneEffect.csv")
    model_path = os.path.join(resolved, "DepMap_Model.csv")

    gene_effect = pd.read_csv(effect_path, index_col=0)
    gene_effect.index.name = "ModelID"
    model = pd.read_csv(model_path)

    # DepMap columns look like "KRAS (3845)"; map the HUGO symbol to the full column name.
    gene_columns = {}
    for column in gene_effect.columns:
        gene_columns[column.split(" (")[0].strip().upper()] = column

    bundle = {
        "gene_effect": gene_effect,
        "model": model,
        "gene_columns": gene_columns,
        "data_lake_path": resolved,
        "provenance": [
            f"CRISPR gene effect: {_file_provenance(effect_path)}; "
            f"{gene_effect.shape[0]} models x {gene_effect.shape[1]} genes (Chronos-corrected)",
            f"Model annotation: {_file_provenance(model_path)}; {model.shape[0]} models",
        ],
    }
    _DEPMAP_CACHE[resolved] = bundle
    return bundle


def _load_expression(data_lake_path: str | None = None) -> dict:
    """Load and cache the DepMap protein-coding expression matrix (log2 TPM+1)."""
    import pandas as pd

    resolved = _resolve_data_lake(data_lake_path)
    if resolved in _EXPRESSION_CACHE:
        return _EXPRESSION_CACHE[resolved]

    path = os.path.join(resolved, "DepMap_OmicsExpressionProteinCodingGenesTPMLogp1.csv")
    if not os.path.exists(path):
        raise FileNotFoundError(f"Expression matrix not found at {path}; confounder analysis needs it.")

    expression = pd.read_csv(path, index_col=0)
    expression.index.name = "ModelID"
    bundle = {
        "expression": expression,
        "gene_columns": {c.split(" (")[0].strip().upper(): c for c in expression.columns},
        "provenance": (
            f"Expression: {_file_provenance(path)}; {expression.shape[0]} models x {expression.shape[1]} genes "
            "(log2 TPM+1)"
        ),
    }
    _EXPRESSION_CACHE[resolved] = bundle
    return bundle


def _select_cancer_models(model_df, cancer_type: str):
    """Subset the DepMap model table to a cancer context.

    Matches ``cancer_type`` as a case-insensitive substring against OncotreeLineage,
    OncotreePrimaryDisease, OncotreeSubtype and OncotreeCode. ``"pan-cancer"`` / ``"all"``
    returns every model.
    """
    if str(cancer_type).strip().lower() in {"pan-cancer", "pancancer", "all", "any"}:
        return model_df.copy(), "pan-cancer (no lineage filter)"

    query = str(cancer_type).strip().lower()
    # "pancreatic cancer" -> "pancrea" so that lineage/disease/subtype spellings all match.
    stem = re.sub(r"\s*(cancer|carcinoma|tumou?r|adenocarcinoma|neoplasm)s?\s*$", "", query).strip()
    stem = stem[:7] if len(stem) > 7 else stem

    columns = ["OncotreeLineage", "OncotreePrimaryDisease", "OncotreeSubtype", "OncotreeCode"]
    mask = None
    for column in columns:
        if column not in model_df.columns:
            continue
        hit = model_df[column].astype(str).str.contains(stem, case=False, na=False, regex=False)
        mask = hit if mask is None else (mask | hit)

    subset = model_df[mask].copy() if mask is not None else model_df.iloc[0:0].copy()
    lineages = sorted(subset["OncotreeLineage"].dropna().unique().tolist()) if len(subset) else []
    return subset, f"matched '{cancer_type}' (stem '{stem}') -> lineages {lineages}"


def _mutation_status_from_local_file(gene: str, data_lake_path: str) -> dict | None:
    """Read mutation status from a local DepMap OmicsSomaticMutations.csv if the user has one."""
    import pandas as pd

    for filename in ("DepMap_OmicsSomaticMutations.csv", "OmicsSomaticMutations.csv"):
        path = os.path.join(data_lake_path, filename)
        if not os.path.exists(path):
            continue
        mutations = pd.read_csv(path, usecols=lambda c: c in {"ModelID", "HugoSymbol", "ProteinChange", "VariantInfo"})
        hits = mutations[mutations["HugoSymbol"].astype(str).str.upper() == gene.upper()]
        variants = {}
        for _, row in hits.iterrows():
            variants.setdefault(row["ModelID"], []).append(str(row.get("ProteinChange", "NA")))
        return {
            "mutant_variants": {k: ",".join(sorted(set(v))) for k, v in variants.items()},
            "profiled_models": set(mutations["ModelID"].unique()),
            "key": "ModelID",
            "source": f"local {filename} ({_file_provenance(path)})",
        }
    return None


def _mutation_status_from_cbioportal(gene: str) -> dict | None:
    """Fetch cell-line mutation calls for ``gene`` from the cBioPortal CCLE study.

    Returns mutant variants and the set of profiled samples keyed by normalised cell line
    name, so that "profiled but not mutated" (= wild type) can be distinguished from
    "never profiled" (= unknown).
    """
    cache_key = gene.upper()
    if cache_key in _MUTATION_CACHE:
        return _MUTATION_CACHE[cache_key]

    gene_info = _http_get(f"{CBIOPORTAL_API}/genes/{urllib.parse.quote(gene.upper())}")
    if not gene_info or "entrezGeneId" not in gene_info:
        return None

    mutations = _http_get(
        f"{CBIOPORTAL_API}/molecular-profiles/{CCLE_STUDY_ID}_mutations/mutations",
        params={
            "sampleListId": f"{CCLE_STUDY_ID}_all",
            "entrezGeneId": gene_info["entrezGeneId"],
            "projection": "SUMMARY",
        },
        timeout=90,
    )
    samples = _http_get(f"{CBIOPORTAL_API}/studies/{CCLE_STUDY_ID}/samples", timeout=90)
    if mutations is None or samples is None:
        return None

    variants: dict[str, list[str]] = {}
    for record in mutations:
        key = _normalize_cell_line_name(str(record.get("sampleId", "")).split("_")[0])
        change = record.get("proteinChange") or record.get("mutationType") or "mutated"
        variants.setdefault(key, []).append(str(change))

    result = {
        "mutant_variants": {k: ",".join(sorted(set(v))) for k, v in variants.items()},
        "profiled_models": {_normalize_cell_line_name(str(s["sampleId"]).split("_")[0]) for s in samples},
        "key": "normalized_name",
        "source": (
            f"cBioPortal {CCLE_STUDY_ID} mutation calls "
            f"({len(mutations)} {gene.upper()} variants across {len(samples)} profiled cell lines)"
        ),
    }
    _MUTATION_CACHE[cache_key] = result
    return result


def _annotate_mutation_status(models, gene: str, data_lake_path: str, mutation_csv_path: str | None = None):
    """Attach MUT / WT / UNKNOWN status for ``gene`` to a DepMap model subset.

    Precedence: user-supplied CSV > local DepMap somatic mutation file > cBioPortal CCLE API.
    """
    import pandas as pd

    models = models.copy()
    models["_norm_name"] = models["StrippedCellLineName"].map(_normalize_cell_line_name)

    annotation = None
    if mutation_csv_path:
        table = pd.read_csv(mutation_csv_path)
        required = {"ModelID", "HugoSymbol"}
        if not required.issubset(table.columns):
            raise ValueError(f"{mutation_csv_path} must contain columns {sorted(required)}")
        hits = table[table["HugoSymbol"].astype(str).str.upper() == gene.upper()]
        annotation = {
            "mutant_variants": {r["ModelID"]: str(r.get("ProteinChange", "mutated")) for _, r in hits.iterrows()},
            "profiled_models": set(table["ModelID"].unique()),
            "key": "ModelID",
            "source": f"user-supplied mutation table {mutation_csv_path}",
        }
    if annotation is None:
        annotation = _mutation_status_from_local_file(gene, data_lake_path)
    if annotation is None:
        annotation = _mutation_status_from_cbioportal(gene)
    if annotation is None:
        raise RuntimeError(
            f"Could not determine {gene} mutation status: no local mutation table and the cBioPortal "
            "API is unreachable. Supply `mutation_csv_path` with columns ModelID,HugoSymbol[,ProteinChange]."
        )

    key_column = "ModelID" if annotation["key"] == "ModelID" else "_norm_name"

    def classify(key):
        if key in annotation["mutant_variants"]:
            return "MUT"
        return "WT" if key in annotation["profiled_models"] else "UNKNOWN"

    models["MutationStatus"] = models[key_column].map(classify)
    models["Variant"] = models[key_column].map(lambda k: annotation["mutant_variants"].get(k, ""))
    return models, annotation["source"]


def _benjamini_hochberg(pvalues):
    """Return BH-FDR adjusted q-values for a 1-D array of p-values."""
    import numpy as np

    pvalues = np.asarray(pvalues, dtype=float)
    n = len(pvalues)
    order = np.argsort(pvalues)
    ranked = pvalues[order]
    qvalues = ranked * n / (np.arange(n) + 1)
    qvalues = np.minimum.accumulate(qvalues[::-1])[::-1]
    out = np.empty(n, dtype=float)
    out[order] = np.clip(qvalues, 0, 1)
    return out


def _parse_gene_list(genes) -> list[str]:
    """Accept a list, a comma separated string, or a newline separated string of gene symbols."""
    if genes is None:
        return []
    if isinstance(genes, str):
        parts = re.split(r"[,\s]+", genes)
    else:
        parts = list(genes)
    return [str(g).strip().upper() for g in parts if str(g).strip()]


# ---------------------------------------------------------------------------
# Tool 1: data-driven synthetic lethality discovery
# ---------------------------------------------------------------------------
def discover_synthetic_lethal_candidates(
    cancer_type: str,
    target_mutation: str,
    data_lake_path: str | None = None,
    mutation_csv_path: str | None = None,
    top_n: int = 20,
    p_threshold: float = 0.05,
    fdr_threshold: float = 0.25,
    min_effect_difference: float = -0.2,
    max_mutant_mean_effect: float = -0.3,
    exclude_pan_essential: bool = True,
    output_csv_path: str | None = None,
) -> str:
    """Discover candidate synthetic-lethal partner genes of a driver mutation in a cancer type.

    Cell lines of the requested cancer type are split into mutant and wild-type groups for
    ``target_mutation``, and every gene in the DepMap CRISPR knockout screen is tested with a
    Welch t-test for stronger dependency (more negative Chronos gene effect) in the mutant group.
    Pan-essential genes, weak dependencies and low-powered comparisons are flagged or removed so
    that the surviving candidates represent genotype-selective vulnerabilities rather than
    generic essential genes.

    Parameters
    ----------
    cancer_type : str
        Cancer context to analyse, e.g. "Pancreatic Cancer", "Lung", "Colorectal Adenocarcinoma".
        Use "pan-cancer" to analyse all lineages together.
    target_mutation : str
        HUGO symbol of the mutated driver gene defining the two groups, e.g. "KRAS", "TP53", "BRCA1".
    data_lake_path : str, optional
        Directory holding DepMap_CRISPRGeneEffect.csv and DepMap_Model.csv
        (default: "./data/biomni_data/data_lake").
    mutation_csv_path : str, optional
        Path to a mutation table (columns: ModelID, HugoSymbol[, ProteinChange]) overriding the
        default mutation source. When omitted, a local DepMap somatic mutation file is used if
        present, otherwise mutation calls are fetched from the cBioPortal CCLE study.
    top_n : int, optional
        Number of top candidates to report in detail (default: 20).
    p_threshold : float, optional
        Uncorrected Welch t-test p-value cutoff (default: 0.05).
    fdr_threshold : float, optional
        Benjamini-Hochberg q-value cutoff (default: 0.25). Set to 1.0 to disable FDR filtering.
    min_effect_difference : float, optional
        Required difference (mutant mean - wild-type mean) in gene effect; must be negative
        (default: -0.2, i.e. the mutant group must be at least 0.2 more depleted).
    max_mutant_mean_effect : float, optional
        The mutant group mean gene effect must be below this value so that statistically
        significant but biologically irrelevant genes are dropped (default: -0.3).
    exclude_pan_essential : bool, optional
        Drop genes depleted in >80% of all screened cell lines (common essentials) (default: True).
    output_csv_path : str, optional
        If given, the full ranked candidate table is written to this CSV path.

    Returns
    -------
    str
        A research log containing the cohort composition, provenance, QC warnings, and a ranked
        candidate table with effect sizes, p/q-values, selectivity and pan-essentiality metrics.

    """
    import numpy as np
    import pandas as pd
    from scipy import stats

    log = [
        "=" * 78,
        f"SYNTHETIC LETHALITY DISCOVERY - {target_mutation.upper()} in {cancer_type}",
        f"Run at {datetime.now(tz=_UTC).strftime('%Y-%m-%d %H:%M UTC')}",
        "=" * 78,
        "",
    ]

    try:
        bundle = _load_depmap(data_lake_path)
    except FileNotFoundError as e:
        return f"FAILURE: {e}"

    gene_effect = bundle["gene_effect"]
    models = bundle["model"]

    # --- Step 1: cohort definition -----------------------------------------------------------
    cohort, match_note = _select_cancer_models(models, cancer_type)
    cohort = cohort[cohort["ModelID"].isin(gene_effect.index)]
    if len(cohort) == 0:
        return (
            f"FAILURE: no DepMap cell line with CRISPR data matched cancer_type='{cancer_type}'. "
            f"Available lineages: {sorted(models['OncotreeLineage'].dropna().unique().tolist())}"
        )

    log.append("STEP 1 | Cohort definition")
    log.append(f"  Cancer context: {match_note}")
    log.append(f"  Cell lines with CRISPR gene-effect data: {len(cohort)}")

    # --- Step 2: mutation stratification -----------------------------------------------------
    try:
        cohort, mutation_source = _annotate_mutation_status(
            cohort, target_mutation, bundle["data_lake_path"], mutation_csv_path
        )
    except (RuntimeError, ValueError) as e:
        return f"FAILURE: {e}"

    mutant_ids = cohort.loc[cohort["MutationStatus"] == "MUT", "ModelID"].tolist()
    wildtype_ids = cohort.loc[cohort["MutationStatus"] == "WT", "ModelID"].tolist()
    unknown_ids = cohort.loc[cohort["MutationStatus"] == "UNKNOWN", "ModelID"].tolist()

    variant_counts = cohort.loc[cohort["MutationStatus"] == "MUT", "Variant"].value_counts().head(6).to_dict()

    log.append("")
    log.append("STEP 2 | Mutation stratification")
    log.append(f"  Mutation source: {mutation_source}")
    log.append(f"  {target_mutation.upper()}-mutant lines : {len(mutant_ids)}")
    log.append(f"  {target_mutation.upper()}-wild-type lines: {len(wildtype_ids)}")
    log.append(f"  Not profiled (excluded)   : {len(unknown_ids)}")
    if variant_counts:
        log.append(f"  Variant spectrum: {variant_counts}")

    if len(mutant_ids) < 3 or len(wildtype_ids) < 3:
        log.append("")
        log.append(
            f"FAILURE: insufficient group sizes (mutant n={len(mutant_ids)}, wild-type n={len(wildtype_ids)}); "
            "at least 3 per group are required. Broaden `cancer_type` (e.g. 'pan-cancer') or supply a "
            "mutation table with wider coverage."
        )
        return "\n".join(log)

    # --- Step 3: differential dependency testing ---------------------------------------------
    mutant_matrix = gene_effect.loc[mutant_ids]
    wildtype_matrix = gene_effect.loc[wildtype_ids]
    usable_genes = mutant_matrix.columns[(mutant_matrix.notna().sum() >= 3) & (wildtype_matrix.notna().sum() >= 3)]
    mutant_matrix = mutant_matrix[usable_genes]
    wildtype_matrix = wildtype_matrix[usable_genes]

    tstat, pvalue = stats.ttest_ind(
        mutant_matrix.values, wildtype_matrix.values, axis=0, equal_var=False, nan_policy="omit"
    )
    mutant_mean = np.asarray(mutant_matrix.mean())
    wildtype_mean = np.asarray(wildtype_matrix.mean())
    pooled_sd = np.sqrt((np.asarray(mutant_matrix.std(ddof=1)) ** 2 + np.asarray(wildtype_matrix.std(ddof=1)) ** 2) / 2)
    with np.errstate(divide="ignore", invalid="ignore"):
        cohens_d = np.where(pooled_sd > 0, (mutant_mean - wildtype_mean) / pooled_sd, np.nan)

    pvalue = np.asarray(pvalue, dtype=float)
    pvalue = np.where(np.isfinite(pvalue), pvalue, 1.0)

    # Pan-essentiality across the whole screen (all lineages), used as a falsification filter.
    depleted_fraction = (gene_effect[usable_genes] < DEPLETION_THRESHOLD).sum() / gene_effect[
        usable_genes
    ].notna().sum()
    # Selectivity inside the cohort: fraction of mutant lines that are actually dependent.
    mutant_dependent_fraction = (mutant_matrix < DEPLETION_THRESHOLD).sum() / mutant_matrix.notna().sum()
    wildtype_dependent_fraction = (wildtype_matrix < DEPLETION_THRESHOLD).sum() / wildtype_matrix.notna().sum()

    results = pd.DataFrame(
        {
            "gene": [c.split(" (")[0] for c in usable_genes],
            "mutant_mean_effect": mutant_mean,
            "wildtype_mean_effect": wildtype_mean,
            "effect_difference": mutant_mean - wildtype_mean,
            "cohens_d": cohens_d,
            "t_statistic": np.asarray(tstat, dtype=float),
            "p_value": pvalue,
            "q_value": _benjamini_hochberg(pvalue),
            "pct_mutant_dependent": np.asarray(mutant_dependent_fraction) * 100,
            "pct_wildtype_dependent": np.asarray(wildtype_dependent_fraction) * 100,
            "pct_all_lines_dependent": np.asarray(depleted_fraction) * 100,
        }
    )
    results["pan_essential"] = results["pct_all_lines_dependent"] >= PAN_ESSENTIAL_FRACTION * 100

    log.append("")
    log.append("STEP 3 | Differential dependency testing (Welch t-test, mutant vs wild-type)")
    log.append(f"  Genes tested: {len(results)}")
    log.append(f"  Nominally significant (p < {p_threshold}): {(results['p_value'] < p_threshold).sum()}")
    log.append(f"  BH q < {fdr_threshold}: {(results['q_value'] < fdr_threshold).sum()}")

    # --- Step 4: filtering / QC ---------------------------------------------------------------
    selected = results[
        (results["p_value"] < p_threshold)
        & (results["q_value"] < fdr_threshold)
        & (results["effect_difference"] <= min_effect_difference)
        & (results["mutant_mean_effect"] <= max_mutant_mean_effect)
    ].copy()
    n_before_pan_essential = len(selected)
    if exclude_pan_essential:
        selected = selected[~selected["pan_essential"]]

    selected = selected.sort_values(["effect_difference", "p_value"]).reset_index(drop=True)
    selected.insert(0, "rank", np.arange(1, len(selected) + 1))

    log.append("")
    log.append("STEP 4 | Filtering and quality control")
    log.append(f"  After effect-size + significance filters: {n_before_pan_essential}")
    log.append(
        f"  Pan-essential genes removed (depleted in >{PAN_ESSENTIAL_FRACTION:.0%} of all "
        f"{gene_effect.shape[0]} lines): {n_before_pan_essential - len(selected)}"
    )
    log.append(f"  Final candidate count: {len(selected)}")

    # Positive control: the driver itself should be a selective dependency in its own mutant lines.
    control_rows = results[results["gene"].str.upper() == target_mutation.upper()]
    if len(control_rows):
        control = control_rows.iloc[0]
        passed = (control["effect_difference"] <= min_effect_difference) and (control["p_value"] < p_threshold)
        in_final = target_mutation.upper() in set(selected["gene"].str.upper())
        status = "PASS" if passed else "FAIL"
        log.append(
            f"  Positive control ({status}): {target_mutation.upper()} itself shows "
            f"delta={control['effect_difference']:.3f}, p={control['p_value']:.2e}, q={control['q_value']:.3f} "
            f"(mutant {control['mutant_mean_effect']:.3f} vs wild type {control['wildtype_mean_effect']:.3f})."
        )
        if passed and not in_final:
            log.append(
                "    Oncogene addiction is reproduced, but the control does not survive the FDR/effect filters - "
                "expected with a small wild-type group, and a reminder that the q-value cutoff is conservative here."
            )
        elif not passed:
            log.append(
                f"    WARNING: {target_mutation.upper()} is not recovered as a selective dependency of its own "
                "mutant lines. The mutation calls or the cohort may be mis-specified; treat all candidates below "
                "as unreliable until this is resolved."
            )
    else:
        log.append(f"  Positive control unavailable: {target_mutation.upper()} is not in the CRISPR library.")

    warnings = []
    if len(wildtype_ids) < 8:
        warnings.append(
            f"Low statistical power: only {len(wildtype_ids)} wild-type lines. p-values are unstable and a "
            "single outlier line can drive a candidate; replicate in a pan-cancer or external cohort."
        )
    if len(unknown_ids):
        warnings.append(
            f"{len(unknown_ids)} lines lacked mutation calls and were excluded rather than assumed wild type."
        )
    lineages = cohort["OncotreeLineage"].dropna().nunique()
    if lineages > 1:
        warnings.append(
            f"The cohort spans {lineages} lineages; lineage composition may confound the mutant/WT contrast."
        )
    warnings.append(
        "Single-perturbation correlative evidence only - no double-knockout or isogenic validation is "
        "included, so these are hypotheses, not confirmed SL interactions."
    )

    log.append("")
    log.append("QC WARNINGS (falsification checklist)")
    for warning in warnings:
        log.append(f"  - {warning}")

    # --- Step 5: report -----------------------------------------------------------------------
    log.append("")
    log.append(f"TOP {min(top_n, len(selected))} CANDIDATE SYNTHETIC LETHAL PARTNERS")
    log.append(
        f"{'rank':<5}{'gene':<12}{'MUT mean':>10}{'WT mean':>10}{'delta':>9}"
        f"{'d':>7}{'p':>11}{'q':>9}{'%MUT dep':>10}{'%all dep':>10}"
    )
    log.append("-" * 93)
    for _, row in selected.head(top_n).iterrows():
        log.append(
            f"{int(row['rank']):<5}{row['gene']:<12}{row['mutant_mean_effect']:>10.3f}"
            f"{row['wildtype_mean_effect']:>10.3f}{row['effect_difference']:>9.3f}"
            f"{row['cohens_d']:>7.2f}{row['p_value']:>11.2e}{row['q_value']:>9.3f}"
            f"{row['pct_mutant_dependent']:>9.0f}%{row['pct_all_lines_dependent']:>9.0f}%"
        )

    candidate_genes = [g for g in selected["gene"].head(top_n).tolist() if g.upper() != target_mutation.upper()]
    log.append("")
    log.append(f"CANDIDATE_GENES: {', '.join(candidate_genes)}")
    log.append(
        "  (pass this list to validate_sl_candidates_with_pubmed and analyze_ppi_network_for_sl "
        "for literature and protein-network evidence)"
    )

    if output_csv_path:
        os.makedirs(os.path.dirname(os.path.abspath(output_csv_path)), exist_ok=True)
        selected.to_csv(output_csv_path, index=False)
        log.append(f"  Full ranked table written to {output_csv_path} ({len(selected)} rows)")

    log.append("")
    log.append("PROVENANCE")
    for item in bundle["provenance"]:
        log.append(f"  - {item}")
    log.append(f"  - Mutation calls: {mutation_source}")
    log.append(
        f"  - Statistics: Welch two-sample t-test on Chronos gene effect, BH-FDR across {len(results)} genes; "
        f"depletion threshold {DEPLETION_THRESHOLD}, pan-essential cutoff {PAN_ESSENTIAL_FRACTION:.0%}"
    )
    log.append(
        f"  - Mutant models: {', '.join(cohort.loc[cohort['MutationStatus'] == 'MUT', 'StrippedCellLineName'][:40])}"
    )
    log.append(
        f"  - Wild-type models: {', '.join(cohort.loc[cohort['MutationStatus'] == 'WT', 'StrippedCellLineName'])}"
    )
    return "\n".join(log)


# ---------------------------------------------------------------------------
# Tool 2: PubMed literature validation
# ---------------------------------------------------------------------------
_SUPPORT_TERMS = [
    "synthetic lethal",
    "synthetic lethality",
    "synthetically lethal",
    "selective dependency",
    "collateral vulnerability",
    "sensitizes",
    "sensitized",
    "sensitivity to",
    "vulnerability",
    "co-dependency",
    "codependency",
    "essential in",
    "required for the survival",
    "growth inhibition",
    "impaired proliferation",
]
_REFUTATION_TERMS = [
    "not synthetic lethal",
    "no synthetic lethal",
    "failed to replicate",
    "failed to reproduce",
    "could not be reproduced",
    "not reproducible",
    "irreproducible",
    "no evidence",
    "did not confirm",
    "was not confirmed",
    "off-target",
    "independent of kras",
    "is dispensable",
    "dispensable for",
    "no significant difference",
    "contrary to",
    "challenge the",
    "call into question",
    "not required for",
    "nonessential",
    "non-essential",
    "is not essential",
    "were not essential",
    "unable to confirm",
    "could not confirm",
    "reevaluation",
    "re-evaluation",
    "reassessment",
    "controversial",
    "questioned",
    "does not depend",
    "no correlation between",
]
_DIRECT_EVIDENCE_TERMS = [
    "crispr screen",
    "shrna screen",
    "rnai screen",
    "isogenic",
    "knockout",
    "knockdown",
    "xenograft",
    "in vivo",
    "combination treatment",
    "double knockout",
]


def _entrez_esearch(query: str, max_papers: int, email: str | None, api_key: str | None) -> list[str]:
    params = {
        "db": "pubmed",
        "term": query,
        "retmax": str(max_papers),
        "retmode": "json",
        "sort": "relevance",
        "tool": "biomni_synthetic_lethality",
    }
    if email:
        params["email"] = email
    if api_key:
        params["api_key"] = api_key
    payload = _http_get(f"{ENTREZ_API}/esearch.fcgi", params=params, timeout=45)
    if not payload:
        return []
    return payload.get("esearchresult", {}).get("idlist", [])


def _entrez_efetch(pmids: list[str], email: str | None, api_key: str | None) -> list[dict]:
    """Fetch PubMed records as XML and parse title / year / journal / abstract."""
    if not pmids:
        return []
    params = {
        "db": "pubmed",
        "id": ",".join(pmids),
        "retmode": "xml",
        "tool": "biomni_synthetic_lethality",
    }
    if email:
        params["email"] = email
    if api_key:
        params["api_key"] = api_key
    xml_text = _http_get(f"{ENTREZ_API}/efetch.fcgi", params=params, timeout=60, as_json=False)
    if not xml_text:
        return []

    try:
        root = ET.fromstring(xml_text)
    except ET.ParseError as e:
        print(f"[synthetic_lethality] failed to parse PubMed XML: {e}")
        return []

    records = []
    for article in root.findall(".//PubmedArticle"):
        pmid_node = article.find(".//PMID")
        title_node = article.find(".//ArticleTitle")

        # Structured abstracts are split over several AbstractText nodes with a Label attribute.
        abstract_parts = []
        for node in article.findall(".//Abstract/AbstractText"):
            text = "".join(node.itertext()).strip()
            if not text:
                continue
            label = node.get("Label")
            abstract_parts.append(f"{label}: {text}" if label else text)

        year = None
        for path in (".//ArticleDate/Year", ".//PubDate/Year", ".//PubMedPubDate/Year"):
            node = article.find(path)
            if node is not None and node.text:
                year = node.text
                break
        if year is None:
            medline_date = article.find(".//PubDate/MedlineDate")
            if medline_date is not None and medline_date.text:
                match = re.search(r"\d{4}", medline_date.text)
                year = match.group(0) if match else None

        journal_node = article.find(".//Journal/Title")
        pub_types = [t.text for t in article.findall(".//PublicationType") if t.text]

        records.append(
            {
                "pmid": pmid_node.text if pmid_node is not None else "",
                "title": "".join(title_node.itertext()).strip() if title_node is not None else "",
                "abstract": " ".join(abstract_parts),
                "year": int(year) if year and year.isdigit() else None,
                "journal": journal_node.text if journal_node is not None else "",
                "publication_types": pub_types,
            }
        )
    return records


def _score_literature(records: list[dict], candidate_gene: str, mutated_gene: str, disease: str) -> dict:
    """Score literature support for one candidate on a 0-100 scale with an explicit breakdown."""
    support_hits, refute_hits, direct_hits = [], [], []
    disease_hits = 0
    cooccurrence = 0
    recent = 0
    current_year = datetime.now(tz=_UTC).year
    disease_tokens = [t for t in re.split(r"[^a-z]+", disease.lower()) if len(t) > 4]

    for record in records:
        text = f"{record['title']} {record['abstract']}".lower()
        if candidate_gene.lower() in text and mutated_gene.lower() in text:
            cooccurrence += 1
        if disease_tokens and any(token[:7] in text for token in disease_tokens):
            disease_hits += 1
        if record["year"] and record["year"] >= current_year - 5:
            recent += 1
        for term in _SUPPORT_TERMS:
            if term in text:
                support_hits.append((record["pmid"], term))
        for term in _REFUTATION_TERMS:
            if term in text:
                refute_hits.append((record["pmid"], term))
        for term in _DIRECT_EVIDENCE_TERMS:
            if term in text:
                direct_hits.append((record["pmid"], term))

    n = len(records)
    volume_score = min(25, n * 5)
    cooccurrence_score = min(25, cooccurrence * 8)
    support_score = min(20, len({p for p, _ in support_hits}) * 7)
    direct_score = min(15, len({p for p, _ in direct_hits}) * 5)
    context_score = min(10, disease_hits * 4)
    recency_score = min(5, recent * 2)
    refutation_penalty = min(30, len({p for p, _ in refute_hits}) * 10)

    total = max(
        0,
        volume_score
        + cooccurrence_score
        + support_score
        + direct_score
        + context_score
        + recency_score
        - refutation_penalty,
    )

    if n == 0:
        interpretation = "UNEXPLORED - no PubMed record links this pair; novel but entirely unvalidated."
    elif refutation_penalty >= 20:
        interpretation = "CONTESTED - the literature contains explicit negative or non-replication statements."
    elif refutation_penalty > 0 and total >= 45:
        interpretation = (
            "SUPPORTED BUT CONTESTED - widely reported, yet at least one paper reports a negative or "
            "non-replication result; the interaction must be treated as unresolved, not established."
        )
    elif refutation_penalty > 0:
        interpretation = (
            "CONTESTED-SPARSE - little supporting literature and at least one negative/non-replication "
            "statement; verify the flagged PMIDs before spending experimental effort."
        )
    elif total >= 65:
        interpretation = "WELL SUPPORTED - repeatedly reported, including direct experimental evidence."
    elif total >= 35:
        interpretation = "EMERGING - some supporting reports, but sparse or indirect."
    else:
        interpretation = "WEAK - mentioned in the literature without specific support for this interaction."

    return {
        "n_papers": n,
        "score": total,
        "interpretation": interpretation,
        "breakdown": {
            "volume": volume_score,
            "gene_pair_cooccurrence": cooccurrence_score,
            "sl_language": support_score,
            "direct_experimental": direct_score,
            "disease_context": context_score,
            "recency": recency_score,
            "refutation_penalty": -refutation_penalty,
        },
        "support_terms": sorted({t for _, t in support_hits}),
        "refutation_terms": sorted({t for _, t in refute_hits}),
        "direct_evidence_terms": sorted({t for _, t in direct_hits}),
        "refuting_pmids": sorted({p for p, _ in refute_hits}),
    }


def validate_sl_candidates_with_pubmed(
    disease: str,
    mutated_gene: str,
    candidate_genes: list[str] | str,
    max_papers_per_gene: int = 8,
    min_year: int | None = None,
    include_abstracts: bool = True,
    abstract_chars: int = 700,
    email: str | None = None,
    api_key: str | None = None,
) -> str:
    """Validate synthetic-lethal candidate genes against the PubMed literature via the NCBI Entrez API.

    For every candidate gene a PubMed query combining the disease, the mutated driver gene and the
    candidate is executed with esearch; the matching records are retrieved with efetch and their
    abstracts are parsed from XML. Each candidate receives a 0-100 literature support score whose
    components are reported explicitly, including a penalty for refuting/non-replication language,
    so that plausible-but-fragile candidates are actively down-ranked.

    Parameters
    ----------
    disease : str
        Disease context used in the query, e.g. "Pancreatic cancer".
    mutated_gene : str
        The mutated driver gene, e.g. "KRAS".
    candidate_genes : list[str] | str
        Candidate SL partner genes; a list, or a comma/whitespace separated string.
    max_papers_per_gene : int, optional
        Maximum number of PubMed records to retrieve per candidate (default: 8).
    min_year : int, optional
        Only consider papers published in or after this year.
    include_abstracts : bool, optional
        Include the (truncated) abstract text of the top paper per candidate (default: True).
    abstract_chars : int, optional
        Number of abstract characters to show when ``include_abstracts`` is True (default: 700).
    email : str, optional
        Contact e-mail passed to Entrez, as recommended by NCBI usage policy.
    api_key : str, optional
        NCBI API key; raises the rate limit from 3 to 10 requests/second.

    Returns
    -------
    str
        A research log with, per candidate, the PubMed query used, the number of hits, the support
        score with its breakdown, refuting evidence, and the key papers (PMID, year, journal, title,
        optional abstract), followed by a ranked summary table.

    """
    genes = _parse_gene_list(candidate_genes)
    if not genes:
        return "FAILURE: no candidate genes supplied."

    api_key = api_key or os.getenv("NCBI_API_KEY")
    email = email or os.getenv("NCBI_EMAIL")
    delay = 0.12 if api_key else 0.4  # NCBI: 10 req/s with a key, 3 req/s without

    log = [
        "=" * 78,
        f"PUBMED LITERATURE VALIDATION - {mutated_gene.upper()} SL partners in {disease}",
        f"Candidates: {', '.join(genes)}",
        f"Run at {datetime.now(tz=_UTC).strftime('%Y-%m-%d %H:%M UTC')}",
        "=" * 78,
    ]

    date_filter = f' AND ("{min_year}"[PDAT] : "3000"[PDAT])' if min_year else ""
    summaries = []

    for gene in genes:
        pair = f'"{mutated_gene}"[Title/Abstract] AND "{gene}"[Title/Abstract]{date_filter}'
        # Tier 1 anchors the disease context; tier 2 deliberately hunts for synthetic-lethality
        # claims AND their refutations outside that context; tier 3 is the unrestricted pair.
        tiers = [
            (
                "disease-specific",
                f'("{disease}"[Title/Abstract] OR "{disease}"[MeSH Terms]) AND {pair}',
            ),
            (
                "SL-focused",
                f"{pair} AND (synthetic lethal*[Title/Abstract] OR dependency[Title/Abstract] OR "
                "vulnerability[Title/Abstract] OR sensitiz*[Title/Abstract] OR screen[Title/Abstract])",
            ),
            (
                "refutation-focused",
                f"{pair} AND (nonessential[Title/Abstract] OR non-essential[Title/Abstract] OR "
                "dispensable[Title/Abstract] OR reproducib*[Title/Abstract] OR off-target[Title/Abstract] OR "
                "controvers*[Title/Abstract] OR reevaluat*[Title/Abstract] OR re-evaluat*[Title/Abstract] OR "
                "reassess*[Title/Abstract] OR fail*[Title/Abstract])",
            ),
        ]

        pmids: list[str] = []
        used_queries = []
        for label, query in tiers:
            hits = _entrez_esearch(query, max_papers_per_gene, email, api_key)
            time.sleep(delay)
            new = [p for p in hits if p not in pmids]
            used_queries.append(f"[{label}] {len(hits)} hit(s), {len(new)} new")
            pmids.extend(new)

        broadened = False
        if not pmids:
            # SL interactions are often first reported outside the disease context of interest.
            pmids = _entrez_esearch(pair, max_papers_per_gene, email, api_key)
            used_queries.append(f"[gene-pair fallback] {len(pmids)} hit(s)")
            broadened = True
            time.sleep(delay)

        pmids = pmids[: max_papers_per_gene * 2]
        records = _entrez_efetch(pmids, email, api_key)
        time.sleep(delay)
        scored = _score_literature(records, gene, mutated_gene, disease)
        scored["gene"] = gene
        scored["query"] = pair
        scored["broadened"] = broadened
        summaries.append(scored)

        log.append("")
        log.append(f"### {gene}")
        log.append(f"  Queries: {'; '.join(used_queries)}")
        log.append(f"  Base pair term: {pair}")
        if broadened:
            log.append("  NOTE: no disease-specific or SL-focused paper found; the unrestricted pair query was used.")
        log.append(f"  Records retrieved: {scored['n_papers']}")
        log.append(f"  Literature support score: {scored['score']}/100 -> {scored['interpretation']}")
        log.append(f"  Score breakdown: {json.dumps(scored['breakdown'])}")
        if scored["support_terms"]:
            log.append(f"  Supporting language: {', '.join(scored['support_terms'][:8])}")
        if scored["direct_evidence_terms"]:
            log.append(f"  Direct experimental evidence terms: {', '.join(scored['direct_evidence_terms'][:8])}")
        if scored["refutation_terms"]:
            log.append(
                f"  REFUTING language: {', '.join(scored['refutation_terms'][:8])} "
                f"(PMIDs {', '.join(scored['refuting_pmids'][:5])})"
            )
        else:
            log.append("  REFUTING language: none detected in the retrieved abstracts.")

        for record in records[:3]:
            log.append(
                f"  - PMID {record['pmid']} ({record['year'] or 'n.d.'}, {record['journal']}): {record['title']}"
            )
            if include_abstracts and record["abstract"]:
                abstract = record["abstract"][:abstract_chars]
                suffix = "..." if len(record["abstract"]) > abstract_chars else ""
                log.append(f"      Abstract: {abstract}{suffix}")

    summaries.sort(key=lambda s: s["score"], reverse=True)
    log.append("")
    log.append("LITERATURE SUPPORT RANKING")
    log.append(f"{'gene':<12}{'papers':>8}{'score':>8}  interpretation")
    log.append("-" * 78)
    for summary in summaries:
        log.append(f"{summary['gene']:<12}{summary['n_papers']:>8}{summary['score']:>8}  {summary['interpretation']}")

    log.append("")
    log.append("PROVENANCE")
    log.append(f"  - NCBI Entrez E-utilities ({ENTREZ_API}), db=pubmed, esearch + efetch (XML abstracts)")
    log.append(f"  - API key in use: {'yes' if api_key else 'no (3 requests/second limit)'}")
    log.append(
        "  - Scores are computed by deterministic keyword rules over retrieved abstracts, not by an LLM; "
        "a low score means 'not documented', not 'disproven'."
    )
    log.append(
        "  - Keyword-based refutation detection can fire on a phrase used in an unrelated context; always open "
        "the flagged PMIDs before acting on a CONTESTED verdict."
    )
    return "\n".join(log)


# ---------------------------------------------------------------------------
# Tool 3: STRING protein-protein interaction analysis
# ---------------------------------------------------------------------------
def analyze_ppi_network_for_sl(
    target_gene: str,
    candidate_genes: list[str] | str,
    species: int = 9606,
    required_score: int = 400,
    analyze_shared_partners: bool = True,
    partner_limit: int = 250,
    run_enrichment: bool = True,
) -> str:
    """Analyse STRING-DB protein interactions between a mutated driver gene and SL candidates.

    Physical/functional proximity is evaluated in two ways: (i) the direct STRING edge between the
    driver and each candidate, broken down by evidence channel (experimental, database, co-expression,
    text-mining), and (ii) the number of shared first-degree interaction partners, which captures
    pathway-level ("hop-2") proximity typical of genuine SL pairs that never touch physically.

    Parameters
    ----------
    target_gene : str
        The mutated driver gene, e.g. "KRAS".
    candidate_genes : list[str] | str
        Candidate SL partner genes; a list, or a comma/whitespace separated string.
    species : int, optional
        NCBI taxonomy identifier (default: 9606, Homo sapiens).
    required_score : int, optional
        Minimum STRING combined score, 0-1000 (default: 400 = medium confidence).
    analyze_shared_partners : bool, optional
        Compute shared first-degree neighbours between the driver and each candidate (default: True).
    partner_limit : int, optional
        Maximum number of interaction partners fetched per gene for the shared-neighbour analysis
        (default: 250).
    run_enrichment : bool, optional
        Run STRING functional enrichment over the driver + candidate set (default: True).

    Returns
    -------
    str
        A research log with an identifier mapping report, a direct-interaction table with per-channel
        evidence scores, shared-partner counts, a functional-proximity classification per candidate,
        and the enriched pathways of the joint network.

    """
    genes = _parse_gene_list(candidate_genes)
    if not genes:
        return "FAILURE: no candidate genes supplied."
    target_gene = target_gene.upper()
    all_genes = [target_gene] + [g for g in genes if g != target_gene]

    log = [
        "=" * 78,
        f"STRING PPI NETWORK ANALYSIS - {target_gene} vs {len(genes)} candidate(s)",
        f"Species {species}, minimum combined score {required_score}",
        f"Run at {datetime.now(tz=_UTC).strftime('%Y-%m-%d %H:%M UTC')}",
        "=" * 78,
        "",
    ]

    # --- Step 1: identifier mapping -----------------------------------------------------------
    mapping = _http_get(
        f"{STRING_API}/json/get_string_ids",
        params={"identifiers": "\r".join(all_genes), "species": species, "limit": 1, "caller_identity": "biomni"},
        timeout=60,
    )
    resolved = {}
    if mapping:
        for item in mapping:
            resolved[str(item.get("queryItem", "")).upper()] = item.get("preferredName", "")
    unresolved = [g for g in all_genes if g not in resolved]

    log.append("STEP 1 | Identifier mapping")
    log.append(f"  Resolved {len(resolved)}/{len(all_genes)} symbols in STRING v12.")
    if unresolved:
        log.append(f"  WARNING: unresolved symbols (excluded): {', '.join(unresolved)}")
    if target_gene not in resolved:
        return "\n".join(log + ["", f"FAILURE: the driver gene {target_gene} could not be mapped in STRING."])

    # --- Step 2: direct interactions ----------------------------------------------------------
    network = _http_get(
        f"{STRING_API}/json/network",
        params={
            "identifiers": "\r".join(all_genes),
            "species": species,
            "required_score": required_score,
            "caller_identity": "biomni",
        },
        timeout=90,
    )
    if network is None:
        return "\n".join(log + ["", "FAILURE: the STRING network endpoint is unreachable."])

    edges = {}
    for edge in network:
        a, b = str(edge.get("preferredName_A", "")).upper(), str(edge.get("preferredName_B", "")).upper()
        edges[frozenset((a, b))] = edge

    log.append("")
    log.append("STEP 2 | Direct interaction with the driver gene")
    log.append(f"{'candidate':<12}{'combined':>10}{'experim.':>10}{'database':>10}{'coexpr.':>10}{'textmine':>10}")
    log.append("-" * 62)

    direct = {}
    for gene in genes:
        edge = edges.get(frozenset((target_gene, resolved.get(gene, gene))))
        if edge:
            direct[gene] = {
                "combined": float(edge.get("score", 0)),
                "experimental": float(edge.get("escore", 0)),
                "database": float(edge.get("dscore", 0)),
                "coexpression": float(edge.get("ascore", 0)),
                "textmining": float(edge.get("tscore", 0)),
            }
            d = direct[gene]
            log.append(
                f"{gene:<12}{d['combined']:>10.3f}{d['experimental']:>10.3f}"
                f"{d['database']:>10.3f}{d['coexpression']:>10.3f}{d['textmining']:>10.3f}"
            )
        else:
            direct[gene] = None
            log.append(f"{gene:<12}  (no direct edge at combined score >= {required_score})")

    # --- Step 3: shared interaction partners ---------------------------------------------------
    shared_counts = {}
    if analyze_shared_partners:

        def partners_of(gene: str) -> set:
            payload = _http_get(
                f"{STRING_API}/json/interaction_partners",
                params={
                    "identifiers": gene,
                    "species": species,
                    "required_score": required_score,
                    "limit": partner_limit,
                    "caller_identity": "biomni",
                },
                timeout=90,
            )
            if not payload:
                return set()
            return {str(p.get("preferredName_B", "")).upper() for p in payload}

        target_partners = partners_of(target_gene)
        log.append("")
        log.append("STEP 3 | Shared first-degree interaction partners (pathway-level proximity)")
        log.append(f"  {target_gene} has {len(target_partners)} partners at score >= {required_score}.")
        for gene in genes:
            if gene not in resolved:
                continue
            candidate_partners = partners_of(gene)
            shared = sorted(target_partners & candidate_partners)
            shared_counts[gene] = shared
            preview = ", ".join(shared[:10]) + ("..." if len(shared) > 10 else "")
            log.append(f"  {gene:<10} partners={len(candidate_partners):<4} shared={len(shared):<4} {preview}")

    # --- Step 4: proximity classification ------------------------------------------------------
    log.append("")
    log.append("STEP 4 | Functional proximity classification")
    classifications = {}
    for gene in genes:
        edge = direct.get(gene)
        shared = shared_counts.get(gene, [])
        if edge and edge["experimental"] >= 0.4:
            verdict = "DIRECT-EXPERIMENTAL: physical/experimental interaction reported with the driver."
        elif edge and edge["combined"] >= 0.7:
            verdict = "DIRECT-HIGH-CONFIDENCE: strong STRING edge, largely from databases/text mining."
        elif edge:
            verdict = "DIRECT-WEAK: an edge exists but only at medium confidence."
        elif len(shared) >= 10:
            verdict = f"PATHWAY-PROXIMAL: no direct edge, but {len(shared)} shared partners (same functional module)."
        elif len(shared) >= 3:
            verdict = f"WEAKLY-PROXIMAL: {len(shared)} shared partners."
        else:
            verdict = (
                "DISTANT: no direct edge and few shared partners - either a parallel-pathway SL "
                "(mechanistically plausible) or a statistical artefact; prioritise mechanistic follow-up."
            )
        classifications[gene] = verdict
        log.append(f"  {gene:<12}{verdict}")

    log.append("")
    log.append(
        "  Interpretation note: absence of a PPI edge does NOT refute synthetic lethality. Classic SL pairs "
        "(e.g. BRCA1-PARP1) act through parallel pathways rather than physical binding. A strong text-mining-only "
        "edge, on the other hand, may merely reflect co-citation bias."
    )

    # --- Step 5: functional enrichment ---------------------------------------------------------
    if run_enrichment:
        enrichment = _http_get(
            f"{STRING_API}/json/enrichment",
            params={"identifiers": "\r".join(all_genes), "species": species, "caller_identity": "biomni"},
            timeout=90,
        )
        log.append("")
        log.append("STEP 5 | Functional enrichment of the driver + candidate set")
        if enrichment is None:
            log.append("  Enrichment endpoint unreachable; skipped.")
        elif len(enrichment) == 0:
            log.append(
                "  No enriched term - the candidates do not share an annotated pathway with the driver. "
                "Expected for parallel-pathway SL, but it also means no pathway-level corroboration is available."
            )
        else:
            filtered = [e for e in enrichment if e.get("category") in {"Process", "KEGG", "RCTM", "WikiPathways"}]
            filtered.sort(key=lambda e: float(e.get("fdr", 1)))
            for item in filtered[:10]:
                log.append(
                    f"  [{item.get('category')}] {item.get('description')} "
                    f"(genes {item.get('number_of_genes')}, FDR {float(item.get('fdr', 1)):.2e})"
                )
            if not filtered:
                log.append("  Enriched terms exist but none in the Process/KEGG/Reactome/WikiPathways categories.")

    log.append("")
    log.append("PROVENANCE")
    log.append(f"  - STRING-DB REST API ({STRING_API}), species={species}, required_score={required_score}")
    endpoints = ["get_string_ids", "network"]
    if analyze_shared_partners:
        endpoints.append("interaction_partners")
    if run_enrichment:
        endpoints.append("enrichment")
    log.append(f"  - Endpoints used: {', '.join(endpoints)}")
    log.append(f"  - Query identifiers: {', '.join(all_genes)}")
    return "\n".join(log)


# ---------------------------------------------------------------------------
# Tool 4: confounder / falsification analysis
# ---------------------------------------------------------------------------
def check_dependency_confounders(
    target_mutation: str,
    candidate_genes: list[str] | str,
    cancer_type: str = "pan-cancer",
    data_lake_path: str | None = None,
    mutation_csv_path: str | None = None,
    top_biomarkers: int = 5,
) -> str:
    """Test whether a candidate dependency is really driven by the mutation or by a confounder.

    A gene can look synthetic-lethal with a driver mutation simply because the two groups differ in
    something else. This tool runs three falsification checks per candidate across the whole DepMap
    panel and reports whether the driver hypothesis survives:

    1. Co-dependency: correlation between the candidate's and the driver's gene-effect profiles.
       Genuine members of the driver's pathway co-vary with it (e.g. RAF1 with KRAS, r about 0.39).
    2. Expression biomarker: the genes whose expression best predicts the candidate's dependency.
       If the top biomarker is the candidate's own paralog rather than anything driver-related, the
       dependency is a paralog-loss effect, not synthetic lethality with the driver (the classic
       VPS4A/VPS4B case).
    3. Lineage effect: how concentrated the dependency is in single lineages, which is the usual
       source of a spurious mutant-versus-wild-type difference.

    Parameters
    ----------
    target_mutation : str
        The driver gene whose mutation defined the candidates, e.g. "KRAS".
    candidate_genes : list[str] | str
        Candidate genes to scrutinise; a list, or a comma/whitespace separated string.
    cancer_type : str, optional
        Restrict the lineage effect summary to this context; correlations always use all lineages
        so that the sample size is adequate (default: "pan-cancer").
    data_lake_path : str, optional
        Directory holding the DepMap files (default: "./data/biomni_data/data_lake").
    mutation_csv_path : str, optional
        Optional mutation table overriding the default mutation source.
    top_biomarkers : int, optional
        Number of top expression biomarkers to report per candidate (default: 5).

    Returns
    -------
    str
        A research log with, per candidate, the driver co-dependency, the strongest expression
        biomarkers of the dependency, the paralog check, the lineage concentration, and a verdict of
        DRIVER-CONSISTENT / CONFOUNDED / UNEXPLAINED.

    """
    import numpy as np
    import pandas as pd
    from scipy import stats

    genes = _parse_gene_list(candidate_genes)
    if not genes:
        return "FAILURE: no candidate genes supplied."

    try:
        bundle = _load_depmap(data_lake_path)
        expression_bundle = _load_expression(data_lake_path)
    except FileNotFoundError as e:
        return f"FAILURE: {e}"

    gene_effect = bundle["gene_effect"]
    effect_columns = bundle["gene_columns"]
    expression = expression_bundle["expression"]
    driver = target_mutation.upper()

    log = [
        "=" * 78,
        f"CONFOUNDER ANALYSIS - are these dependencies really about {driver}?",
        f"Candidates: {', '.join(genes)}",
        f"Run at {datetime.now(tz=_UTC).strftime('%Y-%m-%d %H:%M UTC')}",
        "=" * 78,
        "",
    ]

    if driver not in effect_columns:
        return "\n".join(log + [f"FAILURE: {driver} is not in the CRISPR library, so co-dependency cannot be tested."])

    driver_effect = gene_effect[effect_columns[driver]]

    # Reference scale: how strongly do known members of the driver's own pathway co-vary with it?
    reference = {}
    for pathway_gene in ("RAF1", "SHOC2", "MAPK1", "BRAF", "PTPN11", "SOS1"):
        if pathway_gene in effect_columns and pathway_gene != driver:
            y = gene_effect[effect_columns[pathway_gene]]
            mask = driver_effect.notna() & y.notna()
            if mask.sum() > 30:
                reference[pathway_gene] = stats.pearsonr(driver_effect[mask], y[mask])[0]
    if reference:
        best_reference = max(reference.values())
        log.append("REFERENCE SCALE | co-dependency of known MAPK-pathway genes with " + driver)
        log.append("  " + ", ".join(f"{g} r={r:+.3f}" for g, r in sorted(reference.items(), key=lambda x: -x[1])))
        log.append(f"  Strongest reference r = {best_reference:+.3f}; candidates are judged against this scale.")
    else:
        best_reference = 0.3
        log.append(f"REFERENCE SCALE | unavailable; falling back to r = {best_reference:.2f} as the bar.")

    # Lineage context for the concentration check.
    cohort, _ = _select_cancer_models(bundle["model"], cancer_type)
    cohort_ids = set(cohort["ModelID"]) & set(gene_effect.index)
    lineage_by_model = bundle["model"].set_index("ModelID")["OncotreeLineage"]

    shared_models = gene_effect.index.intersection(expression.index)
    expression_matrix = expression.loc[shared_models]
    expression_values = expression_matrix.to_numpy(dtype=float)
    expression_valid = ~np.isnan(expression_values)

    verdicts = {}
    for gene in genes:
        log.append("")
        log.append(f"### {gene}")
        if gene not in effect_columns:
            log.append("  Not present in the CRISPR library; skipped.")
            continue

        candidate_effect = gene_effect[effect_columns[gene]]

        # --- check 1: co-dependency with the driver ---
        mask = driver_effect.notna() & candidate_effect.notna()
        codependency, codependency_p = stats.pearsonr(driver_effect[mask], candidate_effect[mask])
        log.append(
            f"  1) Co-dependency with {driver}: r={codependency:+.3f} (p={codependency_p:.1e}, n={int(mask.sum())}) "
            f"vs strongest pathway reference r={best_reference:+.3f}"
        )

        # --- check 2: expression biomarkers of the dependency ---
        y = candidate_effect.loc[shared_models].to_numpy(dtype=float)
        y_valid = ~np.isnan(y)
        # Pearson correlation of the dependency against every expressed gene, computed from raw
        # sums so that each column uses only the cell lines where both values are present.
        usable = expression_valid & y_valid[:, None]
        counts = usable.sum(axis=0).astype(float)
        x = np.where(usable, expression_values, 0.0)
        y_masked = np.where(usable, y[:, None], 0.0)
        sum_x = x.sum(axis=0)
        sum_y = y_masked.sum(axis=0)
        sum_xy = (x * y_masked).sum(axis=0)
        sum_xx = (x * x).sum(axis=0)
        sum_yy = (y_masked * y_masked).sum(axis=0)
        numerator = counts * sum_xy - sum_x * sum_y
        denominator = np.sqrt(
            np.clip(counts * sum_xx - sum_x**2, 0, None) * np.clip(counts * sum_yy - sum_y**2, 0, None)
        )
        with np.errstate(divide="ignore", invalid="ignore"):
            correlations = np.where(denominator > 0, numerator / denominator, np.nan)
        correlations = np.where(counts >= 100, correlations, np.nan)

        biomarker_series = pd.Series(correlations, index=[c.split(" (")[0] for c in expression_matrix.columns])
        strongest = biomarker_series.reindex(biomarker_series.abs().sort_values(ascending=False).index)
        strongest = strongest[strongest.index != gene].head(top_biomarkers)
        log.append(
            f"  2) Strongest expression biomarkers of {gene} dependency (positive r = high expression -> less depleted):"
        )
        for biomarker, value in strongest.items():
            log.append(f"       {biomarker:<12} r={value:+.3f}")

        driver_expression_r = biomarker_series.get(driver, float("nan"))
        driver_rank = (
            int((biomarker_series.abs() > abs(driver_expression_r)).sum()) + 1
            if pd.notna(driver_expression_r)
            else None
        )
        if driver_rank:
            log.append(
                f"       {driver} expression itself: r={driver_expression_r:+.3f} (rank {driver_rank} of "
                f"{int(biomarker_series.notna().sum())} genes)"
            )

        # Paralog check: scan the whole gene family explicitly, since a paralog that ranks just
        # outside the top biomarkers (VPS4B for VPS4A) still explains the dependency away.
        family_prefixes = {re.sub(r"(?<=\d)[A-Z]$", "", gene), re.sub(r"\d+$", "", gene)}
        family_prefixes = {p for p in family_prefixes if len(p) >= 3 and p != gene}
        family = [
            symbol
            for symbol in biomarker_series.index
            if symbol != gene and any(symbol.startswith(prefix) for prefix in family_prefixes)
        ]
        paralog_hit = None
        if family:
            family_series = biomarker_series.reindex(sorted(set(family))).dropna()
            family_series = family_series.reindex(family_series.abs().sort_values(ascending=False).index)
            shown = ", ".join(f"{symbol} r={value:+.3f}" for symbol, value in family_series.head(4).items())
            log.append(f"       Paralog family ({'/'.join(sorted(family_prefixes))}*): {shown}")
            if len(family_series) and abs(family_series.iloc[0]) >= 0.25:
                paralog_hit = (family_series.index[0], family_series.iloc[0])
        if paralog_hit:
            paralog_rank = int((biomarker_series.abs() > abs(paralog_hit[1])).sum()) + 1
            log.append(
                f"       PARALOG ALERT: {paralog_hit[0]} expression (r={paralog_hit[1]:+.3f}, rank {paralog_rank} of "
                f"{int(biomarker_series.notna().sum())}) predicts this dependency far better than {driver} "
                f"(rank {driver_rank}) - this looks like a paralog-loss dependency, independent of {driver} status."
            )

        # --- check 3: lineage concentration ---
        depleted = candidate_effect[candidate_effect < DEPLETION_THRESHOLD]
        lineage_counts = lineage_by_model.reindex(depleted.index).value_counts()
        total_depleted = int(lineage_counts.sum())
        if total_depleted:
            top_lineage = lineage_counts.index[0]
            top_share = lineage_counts.iloc[0] / total_depleted * 100
            log.append(
                f"  3) Lineage concentration: {total_depleted} dependent lines overall, "
                f"top lineage {top_lineage} holds {top_share:.0f}%"
            )
        else:
            top_share = 0.0
            log.append("  3) Lineage concentration: no line passes the depletion threshold.")

        in_cohort = [m for m in depleted.index if m in cohort_ids]
        log.append(f"     Within '{cancer_type}': {len(in_cohort)} dependent lines")

        # --- verdict ---
        if paralog_hit:
            verdict = (
                f"CONFOUNDED - the dependency tracks {paralog_hit[0]} expression, not {driver} status. "
                "Stratify by that biomarker before any experiment; as it stands this is not a "
                f"{driver} synthetic lethality."
            )
        elif abs(codependency) >= best_reference * 0.5 and codependency > 0:
            verdict = f"DRIVER-CONSISTENT - co-dependency with {driver} is on the scale of its known pathway members."
        elif pd.notna(driver_expression_r) and driver_rank and driver_rank <= 100:
            verdict = f"DRIVER-PLAUSIBLE - {driver} expression is among the strongest predictors of this dependency."
        elif top_share >= 50:
            verdict = (
                f"CONFOUNDED - {top_share:.0f}% of dependent lines come from one lineage; the mutant/wild-type "
                "contrast may just be lineage composition."
            )
        else:
            verdict = (
                f"UNEXPLAINED - no co-dependency with {driver}, no {driver}-related biomarker and no single "
                "lineage explains it. Could still be a genuine parallel-pathway SL, but the mechanism is open."
            )
        verdicts[gene] = verdict
        log.append(f"  VERDICT: {verdict}")

    log.append("")
    log.append("SUMMARY")
    for gene, verdict in verdicts.items():
        log.append(f"  {gene:<12}{verdict.split(' - ')[0]}")

    log.append("")
    log.append("PROVENANCE")
    for item in bundle["provenance"]:
        log.append(f"  - {item}")
    log.append(f"  - {expression_bundle['provenance']}")
    log.append(
        "  - Correlations are Pearson across all DepMap lines with both measurements; a confounded verdict is a "
        "reason to stratify the analysis, not proof that the interaction is absent."
    )
    return "\n".join(log)


# ---------------------------------------------------------------------------
# Tool 5: integrated evidence dossier
# ---------------------------------------------------------------------------
def _parse_candidate_line(discovery_report: str) -> list[str]:
    for line in discovery_report.splitlines():
        if line.startswith("CANDIDATE_GENES:"):
            return _parse_gene_list(line.split(":", 1)[1])
    return []


def _parse_discovery_table(discovery_report: str) -> dict:
    """Recover the numeric candidate table from a discovery report for integrated scoring."""
    stats = {}
    for line in discovery_report.splitlines():
        match = re.match(
            r"^\s*(\d+)\s+([A-Z0-9\-\.]+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+"
            r"(-?\d+\.\d+)\s+([\d.eE+-]+)\s+(\d+\.\d+)\s+(\d+)%\s+(\d+)%\s*$",
            line,
        )
        if match:
            stats[match.group(2)] = {
                "rank": int(match.group(1)),
                "mutant_mean_effect": float(match.group(3)),
                "wildtype_mean_effect": float(match.group(4)),
                "effect_difference": float(match.group(5)),
                "cohens_d": float(match.group(6)),
                "p_value": float(match.group(7)),
                "q_value": float(match.group(8)),
                "pct_mutant_dependent": float(match.group(9)),
                "pct_all_lines_dependent": float(match.group(10)),
            }
    return stats


def _parse_literature_scores(literature_report: str) -> dict:
    scores = {}
    current = None
    for line in literature_report.splitlines():
        header = re.match(r"^### ([A-Z0-9\-\.]+)\s*$", line)
        if header:
            current = header.group(1)
            continue
        if current:
            score_match = re.search(r"Literature support score: (\d+)/100 -> (.+)$", line)
            if score_match:
                scores[current] = {
                    "score": int(score_match.group(1)),
                    "interpretation": score_match.group(2).strip(),
                    "refuted": False,
                }
            if line.strip().startswith("REFUTING language:") and "none detected" not in line:
                scores.setdefault(current, {"score": 0, "interpretation": "", "refuted": False})
                scores[current]["refuted"] = True
    return scores


def _parse_ppi_classification(ppi_report: str) -> dict:
    verdicts = {}
    in_section = False
    for line in ppi_report.splitlines():
        if line.startswith("STEP 4 |"):
            in_section = True
            continue
        if in_section:
            if line.startswith("  Interpretation note") or line.startswith("STEP 5"):
                in_section = False
                continue
            match = re.match(r"^\s{2}([A-Z0-9\-\.]+)\s{2,}([A-Z\-]+):", line)
            if match:
                verdicts[match.group(1)] = match.group(2)
    return verdicts


def generate_sl_evidence_dossier(
    cancer_type: str,
    target_mutation: str,
    top_n: int = 5,
    data_lake_path: str | None = None,
    mutation_csv_path: str | None = None,
    max_papers_per_gene: int = 8,
    required_score: int = 400,
    email: str | None = None,
    api_key: str | None = None,
) -> str:
    """Run the full synthetic lethality pipeline and produce an integrated evidence dossier.

    Chains ``discover_synthetic_lethal_candidates`` (statistical evidence),
    ``validate_sl_candidates_with_pubmed`` (literature evidence) and ``analyze_ppi_network_for_sl``
    (protein-network evidence), then integrates the three evidence streams into a per-candidate
    dossier with a confidence grade, supporting and contradicting evidence, remaining uncertainties,
    a minimal validation experiment and a Go / Hold / No-go recommendation.

    Because the three streams are not independent (all three ultimately derive from partially
    overlapping cell-line and literature resources), evidence is weighted by directness rather than
    simply summed: genotype-selective CRISPR effect size dominates, literature and PPI evidence
    modulate the confidence grade, and explicit refuting evidence can veto a Go recommendation.

    Parameters
    ----------
    cancer_type : str
        Cancer context, e.g. "Pancreatic Cancer".
    target_mutation : str
        Mutated driver gene defining the genotype contrast, e.g. "KRAS".
    top_n : int, optional
        Number of top statistical candidates carried into literature and PPI validation (default: 5).
    data_lake_path : str, optional
        Directory holding the DepMap files (default: "./data/biomni_data/data_lake").
    mutation_csv_path : str, optional
        Optional mutation table overriding the default mutation source.
    max_papers_per_gene : int, optional
        Maximum PubMed records per candidate (default: 8).
    required_score : int, optional
        Minimum STRING combined score (default: 400).
    email : str, optional
        Contact e-mail for the NCBI Entrez API.
    api_key : str, optional
        NCBI API key for a higher Entrez rate limit.

    Returns
    -------
    str
        The full research log of all three stages followed by an integrated evidence dossier
        (confidence grade, supporting/contradicting evidence, proposed validation experiment and
        Go/Hold/No-go recommendation per candidate).

    """
    sections = []

    discovery = discover_synthetic_lethal_candidates(
        cancer_type=cancer_type,
        target_mutation=target_mutation,
        data_lake_path=data_lake_path,
        mutation_csv_path=mutation_csv_path,
        top_n=max(top_n, 10),
    )
    sections.append(discovery)
    if discovery.startswith("FAILURE") or "CANDIDATE_GENES:" not in discovery:
        return "\n\n".join(sections + ["Pipeline stopped: candidate discovery produced no usable candidates."])

    candidates = _parse_candidate_line(discovery)[:top_n]
    discovery_stats = _parse_discovery_table(discovery)

    literature = validate_sl_candidates_with_pubmed(
        disease=cancer_type,
        mutated_gene=target_mutation,
        candidate_genes=candidates,
        max_papers_per_gene=max_papers_per_gene,
        email=email,
        api_key=api_key,
    )
    sections.append(literature)
    literature_scores = _parse_literature_scores(literature)

    ppi = analyze_ppi_network_for_sl(
        target_gene=target_mutation,
        candidate_genes=candidates,
        required_score=required_score,
    )
    sections.append(ppi)
    ppi_verdicts = _parse_ppi_classification(ppi)

    # --- integration ---------------------------------------------------------------------------
    dossier = [
        "=" * 78,
        f"EVIDENCE DOSSIER - {target_mutation.upper()}-driven synthetic lethality in {cancer_type}",
        f"Generated {datetime.now(tz=_UTC).strftime('%Y-%m-%d %H:%M UTC')}",
        "=" * 78,
    ]

    ranked = []
    for gene in candidates:
        stats_row = discovery_stats.get(gene, {})
        literature_row = literature_scores.get(gene, {"score": 0, "interpretation": "not evaluated", "refuted": False})
        ppi_verdict = ppi_verdicts.get(gene, "UNKNOWN")

        effect = stats_row.get("effect_difference", 0.0)
        selectivity = stats_row.get("pct_mutant_dependent", 0.0) - 0.0
        pan_essentiality = stats_row.get("pct_all_lines_dependent", 0.0)

        # Deterministic integrated score: directness-weighted, not a rank average.
        statistical_component = min(45, abs(effect) * 60)
        selectivity_component = min(20, selectivity * 0.2)
        # Novelty must not be punished like refutation: an undocumented pair starts from a neutral
        # baseline, and only explicitly refuting literature subtracts (via the penalty below).
        literature_component = min(25, 8 + literature_row["score"] * 0.2)
        ppi_component = {
            "DIRECT-EXPERIMENTAL": 10,
            "DIRECT-HIGH-CONFIDENCE": 7,
            "DIRECT-WEAK": 4,
            "PATHWAY-PROXIMAL": 6,
            "WEAKLY-PROXIMAL": 3,
            "DISTANT": 0,
        }.get(ppi_verdict, 0)
        penalty = 0
        contradictions = []
        if pan_essentiality >= 70:
            penalty += 25
            contradictions.append(
                f"depleted in {pan_essentiality:.0f}% of all screened lines - largely a common essential gene; "
                "the genotype effect is a modulation of a general dependency and the therapeutic window is doubtful"
            )
        elif pan_essentiality >= 50:
            penalty += 15
            contradictions.append(
                f"depleted in {pan_essentiality:.0f}% of all screened lines - partly a common essential gene, "
                "so the therapeutic window in normal tissue is uncertain"
            )
        if literature_row.get("refuted"):
            penalty += 15
            contradictions.append("published refuting / non-replication statements were detected")
        if stats_row.get("q_value", 1) > 0.1:
            penalty += 5
            contradictions.append(
                f"FDR q={stats_row.get('q_value', float('nan')):.3f} - only weak multiple-testing support"
            )
        if ppi_verdict == "DISTANT":
            contradictions.append("no protein-network proximity to the driver; mechanism unexplained")

        total = max(0, statistical_component + selectivity_component + literature_component + ppi_component - penalty)

        if total >= 70 and not literature_row.get("refuted"):
            grade, recommendation = "A", "GO"
        elif total >= 50:
            grade, recommendation = "B", "GO (with the confirmatory experiment below)"
        elif total >= 30:
            grade, recommendation = "C", "HOLD"
        else:
            grade, recommendation = "D", "NO-GO"

        ranked.append(
            {
                "gene": gene,
                "total": total,
                "grade": grade,
                "recommendation": recommendation,
                "stats": stats_row,
                "literature": literature_row,
                "ppi": ppi_verdict,
                "contradictions": contradictions,
                "components": {
                    "statistical": round(statistical_component, 1),
                    "selectivity": round(selectivity_component, 1),
                    "literature": round(literature_component, 1),
                    "ppi": ppi_component,
                    "penalty": -penalty,
                },
            }
        )

    ranked.sort(key=lambda r: r["total"], reverse=True)

    dossier.append("")
    dossier.append(f"{'rank':<6}{'gene':<12}{'score':>7}{'grade':>7}  recommendation")
    dossier.append("-" * 78)
    for index, item in enumerate(ranked, start=1):
        dossier.append(f"{index:<6}{item['gene']:<12}{item['total']:>7.1f}{item['grade']:>7}  {item['recommendation']}")

    for index, item in enumerate(ranked, start=1):
        stats_row = item["stats"]
        gene = item["gene"]
        dossier.append("")
        dossier.append("-" * 78)
        dossier.append(f"CANDIDATE {index}: {target_mutation.upper()} - {gene}  ({cancer_type})")
        dossier.append("-" * 78)
        dossier.append(
            f"  Conclusion: {gene} is a candidate synthetic-lethal partner of mutant {target_mutation.upper()} "
            f"in {cancer_type} with confidence grade {item['grade']} "
            f"(integrated score {item['total']:.1f}/100)."
        )
        dossier.append(f"  Score composition: {json.dumps(item['components'])}")
        dossier.append("")
        dossier.append("  SUPPORTING EVIDENCE")
        if stats_row:
            dossier.append(
                f"    [DepMap CRISPR] gene effect {stats_row['mutant_mean_effect']:.3f} (mutant) vs "
                f"{stats_row['wildtype_mean_effect']:.3f} (wild type); delta={stats_row['effect_difference']:.3f}, "
                f"Cohen's d={stats_row['cohens_d']:.2f}, p={stats_row['p_value']:.2e}, q={stats_row['q_value']:.3f}; "
                f"{stats_row['pct_mutant_dependent']:.0f}% of mutant lines are dependent."
            )
        dossier.append(
            f"    [PubMed] support score {item['literature']['score']}/100 - {item['literature']['interpretation']}"
        )
        dossier.append(f"    [STRING PPI] {item['ppi']}")
        dossier.append("")
        dossier.append("  CONTRADICTING EVIDENCE / FAILURE MODES")
        if item["contradictions"]:
            for contradiction in item["contradictions"]:
                dossier.append(f"    - {contradiction}")
        else:
            dossier.append("    - No explicit contradiction detected; the main risk remains the small wild-type group.")
        dossier.append("")
        dossier.append("  MINIMAL VALIDATION EXPERIMENT")
        dossier.append(
            f"    1. Isogenic pair: knock out {gene} (2 independent sgRNAs) in a {target_mutation.upper()}-mutant "
            "line and its wild-type counterpart; read out 10-day viability plus caspase-3/7."
        )
        dossier.append(
            f"    2. Rescue: re-express sgRNA-resistant {gene} cDNA; loss of the differential effect confirms "
            "on-target action."
        )
        dossier.append(
            f"    3. Genotype panel: repeat in >=3 mutant and >=3 wild-type lines and test the "
            f"genotype x {gene}-knockout interaction term (two-way ANOVA)."
        )
        dossier.append(
            "    4. Therapeutic window: repeat in a non-transformed line (e.g. HPNE/HPDE) - a comparable "
            "dependency there is a no-go signal."
        )
        dossier.append("    Decision criterion: >=2-fold viability difference (mutant vs wild type), p<0.05, rescued.")
        dossier.append("")
        dossier.append(f"  RECOMMENDATION: {item['recommendation']}")

    dossier.append("")
    dossier.append("UNRESOLVED QUESTIONS / DATA GAPS")
    dossier.append(
        "  - All statistical evidence derives from a single DepMap release; results were not replicated in an "
        "independent screen (e.g. Sanger Project Score) and must not be counted as independent confirmation."
    )
    dossier.append(
        "  - Single-gene knockout dependency is a proxy for synthetic lethality; no combinatorial CRISPR data is used."
    )
    dossier.append(
        "  - Literature scores measure documentation, not truth: an unexplored pair scores low yet may still be real."
    )
    dossier.append(
        "  - Cell-line dependency does not establish a therapeutic window in patients; normal-tissue toxicity is untested."
    )

    sections.append("\n".join(dossier))
    return "\n\n".join(sections)
