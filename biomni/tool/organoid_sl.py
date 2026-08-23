"""Organoid-aware synthetic lethality tools for Biomni.

Cell-line CRISPR screens are the only large-scale dependency resource that exists, but they
systematically miss a class of synthetic lethal interactions. This module makes those blind spots
computable, so that a candidate can be routed to patient-derived organoid (PDO) validation for a
stated reason rather than as a vague "let's also try organoids".

Three blind spots are addressed, each with data that is actually available locally:

1. ``discover_subtype_masked_sl_candidates`` - a pooled mutant-vs-wild-type test averages over
   transcriptional subtypes (PDAC classical vs basal-like). An interaction restricted to one subtype
   cancels out in the pooled analysis and never reaches the candidate list. This tool runs the test
   inside each subtype and reports what the pooled analysis hid.
2. ``assess_organoid_transferability`` - organoid medium supplies WNT/RSPO1, EGF, FGF10, Noggin, a
   TGF-beta receptor inhibitor and a ROCK inhibitor, and the cells grow anchored in Matrigel rather
   than on plastic. Each of those changes which genes are limiting. This tool predicts, per
   candidate, whether the dependency will be masked, enhanced or preserved in an organoid, and says
   how to modify the medium so the experiment can still answer the question.
3. ``design_organoid_sl_experiment`` - turns the above into a concrete PDO protocol with a matched
   normal-organoid arm, which is the therapeutic-window control that no cell-line panel provides.

Honest limitation, stated up front and repeated in every tool's output: the local DepMap snapshot
contains 5 pancreatic organoid models (PANFR series) but **no CRISPR or expression data for them**.
Nothing here measures an organoid dependency. These tools predict where cell-line evidence is
untrustworthy and design the organoid experiment that would settle it.
"""

import os
from datetime import datetime

from biomni.tool.synthetic_lethality import (
    _UTC,
    DEPLETION_THRESHOLD,
    _annotate_mutation_status,
    _benjamini_hochberg,
    _load_depmap,
    _load_expression,
    _parse_gene_list,
    _resolve_data_lake,
    _select_cancer_models,
)

_GTEX_CACHE: dict = {}
_GENETIC_INTERACTION_CACHE: dict = {}

# Representative transcriptional subtype markers. PDAC classical/basal-like follows Moffitt et al.
# (Nat Genet 2015) tumour-specific subtypes; these are widely used marker genes, not the full
# published signature, and the tool reports them so the choice stays auditable.
SUBTYPE_MARKERS = {
    "pancreatic": {
        "classical": [
            "GATA6",
            "TFF1",
            "TFF2",
            "TFF3",
            "AGR2",
            "LGALS4",
            "CEACAM5",
            "CEACAM6",
            "ANXA10",
            "CTSE",
            "MYO1A",
            "ST6GALNAC1",
            "REG4",
            "SPINK4",
            "CLRN3",
            "VSIG2",
            "FAM3D",
        ],
        "basal": [
            "KRT5",
            "KRT6A",
            "KRT14",
            "KRT17",
            "S100A2",
            "SPRR3",
            "SPRR1B",
            "DHRS9",
            "CST6",
            "LY6D",
            "FAM83A",
            "TNS4",
            "GPR87",
            "SERPINB3",
            "SERPINB4",
            "VGLL1",
            "AREG",
        ],
    },
    "colorectal": {
        "classical": ["CDX2", "LGR5", "ASCL2", "AXIN2", "EPHB2", "OLFM4"],
        "basal": ["ZEB1", "VIM", "SNAI2", "TWIST1", "FN1"],
    },
}

# Standard PDAC organoid medium (Boj et al., Cell 2015; Tuveson-lab human complete feeding medium).
# A niche factor supplied in the medium relieves the cell of having to make or activate it, so
# dependencies UPSTREAM of that factor disappear in an organoid while downstream ones survive.
NICHE_FACTORS = {
    "WNT3A + RSPO1": {
        "masked_upstream": ["PORCN", "WLS", "WNT7B", "WNT5A", "WNT3", "WNT10A", "RNF43", "ZNRF3"],
        "preserved_downstream": ["CTNNB1", "TCF7L2", "LEF1", "TCF7", "AXIN1", "APC", "GSK3B", "CSNK1A1"],
        "note": "exogenous WNT3A/RSPO1 rescues loss of WNT ligand production and of RNF43/ZNRF3 turnover",
    },
    "EGF": {
        "masked_upstream": ["EGF", "TGFA", "AREG", "EREG", "HBEGF", "ADAM17"],
        "preserved_downstream": ["EGFR", "ERBB2", "ERBB3", "GRB2", "SOS1", "SHC1"],
        "note": "saturating EGF removes autocrine ligand dependency but the receptor is still required",
    },
    "FGF10": {
        "masked_upstream": ["FGF10", "FGF7", "FGF2", "FGF1"],
        "preserved_downstream": ["FGFR1", "FGFR2", "FGFR3", "FRS2"],
        "note": "exogenous FGF10 substitutes for autocrine FGF ligands",
    },
    "Noggin (BMP inhibitor)": {
        "masked_upstream": ["BMPR1A", "BMPR2", "SMAD1", "SMAD5", "SMAD9", "ID1", "ID2"],
        "preserved_downstream": [],
        "note": "BMP signalling is pharmacologically blocked, so BMP-arm knockouts add little",
    },
    "A83-01 (TGF-beta receptor inhibitor)": {
        "masked_upstream": ["TGFBR1", "TGFBR2", "SMAD2", "SMAD3", "SMAD4", "TGFB1"],
        "preserved_downstream": [],
        "note": "TGF-beta signalling is already inhibited; TGFBR/SMAD knockouts are largely redundant",
    },
    "Y-27632 (ROCK inhibitor)": {
        "masked_upstream": ["ROCK1", "ROCK2"],
        "preserved_downstream": ["RHOA", "MYH9"],
        "note": "ROCK is pharmacologically inhibited during seeding; ROCK knockouts are masked",
    },
    "Matrigel / laminin-rich ECM": {
        "masked_upstream": ["LAMA1", "LAMB1", "LAMC1", "COL4A1", "FN1"],
        "preserved_downstream": ["ITGB1", "ITGA6", "PTK2", "ILK", "YAP1", "TEAD1", "RHOA"],
        "note": "matrix ligands are supplied, but adhesion signalling is engaged more, not less",
    },
}


def _load_gtex(data_lake_path: str | None = None):
    """Load and cache the GTEx tissue median expression table (long format)."""
    import pandas as pd

    resolved = _resolve_data_lake(data_lake_path)
    if resolved in _GTEX_CACHE:
        return _GTEX_CACHE[resolved]
    path = os.path.join(resolved, "gtex_tissue_gene_tpm.parquet")
    if not os.path.exists(path):
        return None
    table = pd.read_parquet(path)
    _GTEX_CACHE[resolved] = table
    return table


def _load_human_genetic_interactions(data_lake_path: str | None = None):
    """Load and cache human (taxid 9606) genetic interactions with HUGO symbols attached."""
    import pandas as pd

    resolved = _resolve_data_lake(data_lake_path)
    if resolved in _GENETIC_INTERACTION_CACHE:
        return _GENETIC_INTERACTION_CACHE[resolved]

    gi_path = os.path.join(resolved, "genetic_interaction.parquet")
    info_path = os.path.join(resolved, "gene_info.parquet")
    if not (os.path.exists(gi_path) and os.path.exists(info_path)):
        return None

    interactions = pd.read_parquet(gi_path)
    interactions = interactions[(interactions.organism_id_a == 9606) & (interactions.organism_id_b == 9606)].copy()
    info = pd.read_parquet(info_path, columns=["gene_id", "gene_name"]).dropna().drop_duplicates("gene_id")
    mapping = dict(zip(info.gene_id, info.gene_name, strict=False))
    interactions["symbol_a"] = interactions.gene_a_id.map(mapping)
    interactions["symbol_b"] = interactions.gene_b_id.map(mapping)
    interactions = interactions.dropna(subset=["symbol_a", "symbol_b"])
    _GENETIC_INTERACTION_CACHE[resolved] = interactions
    return interactions


def _subtype_scores(expression, expression_columns, model_ids, marker_set: dict):
    """Return a per-model (classical - basal) z-score using the given marker sets."""
    import numpy as np

    usable = [m for m in model_ids if m in expression.index]
    if len(usable) < 6:
        return None, []
    block = expression.loc[usable]
    z = (block - block.mean()) / block.std().replace(0, np.nan)

    used = {}
    scores = {}
    for label in ("classical", "basal"):
        columns = [expression_columns[g] for g in marker_set[label] if g in expression_columns]
        used[label] = [g for g in marker_set[label] if g in expression_columns]
        scores[label] = z[columns].mean(axis=1)
    return (scores["classical"] - scores["basal"]), used


def discover_subtype_masked_sl_candidates(
    cancer_type: str,
    target_mutation: str,
    subtype_marker_set: str = "pancreatic",
    data_lake_path: str | None = None,
    mutation_csv_path: str | None = None,
    p_threshold: float = 0.05,
    min_effect_difference: float = -0.2,
    max_mutant_mean_effect: float = -0.3,
    top_n: int = 15,
) -> str:
    """Find synthetic lethal candidates that a pooled cell-line analysis hides by averaging subtypes.

    A mutant-versus-wild-type test run over a whole cell-line panel assumes the interaction is
    uniform across the panel. Carcinomas are not uniform: PDAC splits into classical and basal-like
    transcriptional subtypes with different differentiation programmes, and an interaction present
    only in one subtype is diluted to non-significance when both are pooled. Patient-derived
    organoids retain the subtype of the tumour they came from, so subtype-restricted interactions are
    exactly the ones worth taking to an organoid.

    This tool scores every cell line for subtype from expression, splits the panel, repeats the
    mutation-stratified Welch t-test inside each subtype, and reports genes that are significant in
    one subtype but not in the pooled test.

    Parameters
    ----------
    cancer_type : str
        Cancer context, e.g. "Pancreatic Cancer".
    target_mutation : str
        Driver gene defining the genotype contrast, e.g. "KRAS".
    subtype_marker_set : str, optional
        Which built-in marker set to use: "pancreatic" (classical vs basal-like) or "colorectal"
        (default: "pancreatic").
    data_lake_path : str, optional
        Directory holding the DepMap files (default: "./data/biomni_data/data_lake").
    mutation_csv_path : str, optional
        Optional mutation table overriding the default mutation source.
    p_threshold : float, optional
        Uncorrected Welch p-value cutoff used inside each subtype (default: 0.05).
    min_effect_difference : float, optional
        Required mutant-minus-wild-type gene-effect difference; must be negative (default: -0.2).
    max_mutant_mean_effect : float, optional
        The mutant group mean gene effect must be below this value (default: -0.3).
    top_n : int, optional
        Number of subtype-masked candidates to report (default: 15).

    Returns
    -------
    str
        A research log with the subtype composition of the panel, the per-subtype group sizes, and a
        ranked table of genes significant within a subtype but missed by the pooled analysis, with an
        explicit statement of how weak the evidence is and why an organoid is the right next step.

    """
    import numpy as np
    import pandas as pd
    from scipy import stats

    log = [
        "=" * 78,
        f"SUBTYPE-MASKED SL DISCOVERY - {target_mutation.upper()} in {cancer_type}",
        f"Run at {datetime.now(tz=_UTC).strftime('%Y-%m-%d %H:%M UTC')}",
        "=" * 78,
        "",
        "Premise: a pooled mutant-vs-wild-type test averages over transcriptional subtypes and",
        "cancels interactions that exist in only one of them. Those are the candidates that a",
        "cell-line panel structurally cannot surface, and that a subtype-matched organoid can.",
        "",
    ]

    if subtype_marker_set not in SUBTYPE_MARKERS:
        return f"FAILURE: unknown subtype_marker_set '{subtype_marker_set}'. Available: {list(SUBTYPE_MARKERS)}"

    try:
        bundle = _load_depmap(data_lake_path)
        expression_bundle = _load_expression(data_lake_path)
    except FileNotFoundError as e:
        return f"FAILURE: {e}"

    gene_effect = bundle["gene_effect"]
    cohort, match_note = _select_cancer_models(bundle["model"], cancer_type)
    cohort = cohort[cohort["ModelID"].isin(gene_effect.index)]
    if len(cohort) < 10:
        return f"FAILURE: only {len(cohort)} screened cell lines matched '{cancer_type}'; subtype splitting needs more."

    try:
        cohort, mutation_source = _annotate_mutation_status(
            cohort, target_mutation, bundle["data_lake_path"], mutation_csv_path
        )
    except (RuntimeError, ValueError) as e:
        return f"FAILURE: {e}"

    expression = expression_bundle["expression"]
    expression_columns = expression_bundle["gene_columns"]
    scores, markers_used = _subtype_scores(
        expression, expression_columns, cohort["ModelID"].tolist(), SUBTYPE_MARKERS[subtype_marker_set]
    )
    if scores is None:
        return "FAILURE: too few cell lines have expression data to compute subtype scores."

    cohort = cohort[cohort["ModelID"].isin(scores.index)].copy()
    cohort["subtype_score"] = cohort["ModelID"].map(scores)
    cohort["subtype"] = np.where(cohort["subtype_score"] > 0, "classical-like", "basal-like")

    log.append("STEP 1 | Panel composition")
    log.append(f"  Context: {match_note}")
    log.append(f"  Mutation source: {mutation_source}")
    log.append(
        f"  Subtype markers used: {len(markers_used['classical'])} classical / {len(markers_used['basal'])} basal "
        f"({subtype_marker_set}, Moffitt-style representative markers)"
    )
    composition = cohort.groupby(["subtype", "MutationStatus"]).size().unstack(fill_value=0)
    for subtype, row in composition.iterrows():
        log.append(
            f"  {subtype:<15} mutant={int(row.get('MUT', 0)):<3} wild-type={int(row.get('WT', 0)):<3} "
            f"unknown={int(row.get('UNKNOWN', 0))}"
        )
    log.append("  NOTE: subtype scores are z-scored within this panel, so they rank lines relative to each")
    log.append("  other. They do NOT establish how the panel compares with patient tumours - no patient")
    log.append("  reference cohort is present in the data lake.")

    def run_test(model_ids):
        mutant = [m for m in model_ids if cohort.set_index("ModelID").loc[m, "MutationStatus"] == "MUT"]
        wildtype = [m for m in model_ids if cohort.set_index("ModelID").loc[m, "MutationStatus"] == "WT"]
        if len(mutant) < 3 or len(wildtype) < 2:
            return None, len(mutant), len(wildtype)
        a = gene_effect.loc[mutant]
        b = gene_effect.loc[wildtype]
        usable = a.columns[(a.notna().sum() >= 3) & (b.notna().sum() >= 2)]
        a, b = a[usable], b[usable]
        tstat, pvalue = stats.ttest_ind(a.values, b.values, axis=0, equal_var=False, nan_policy="omit")
        pvalue = np.where(np.isfinite(np.asarray(pvalue, dtype=float)), np.asarray(pvalue, dtype=float), 1.0)
        mutant_mean = np.asarray(a.mean())
        wildtype_mean = np.asarray(b.mean())
        frame = pd.DataFrame(
            {
                "gene": [c.split(" (")[0] for c in usable],
                "mutant_mean": mutant_mean,
                "wildtype_mean": wildtype_mean,
                "difference": mutant_mean - wildtype_mean,
                "p_value": pvalue,
                "q_value": _benjamini_hochberg(pvalue),
            }
        ).set_index("gene")
        return frame, len(mutant), len(wildtype)

    pooled, pooled_mut, pooled_wt = run_test(cohort["ModelID"].tolist())
    if pooled is None:
        return "\n".join(log + ["", "FAILURE: the pooled comparison has too few lines per genotype group."])

    log.append("")
    log.append("STEP 2 | Pooled versus per-subtype testing")
    log.append(f"  Pooled          : mutant n={pooled_mut}, wild-type n={pooled_wt}")

    pan_essential = (gene_effect < DEPLETION_THRESHOLD).sum() / gene_effect.notna().sum()
    pan_essential.index = [c.split(" (")[0] for c in gene_effect.columns]

    findings = []
    for subtype in ("classical-like", "basal-like"):
        ids = cohort.loc[cohort["subtype"] == subtype, "ModelID"].tolist()
        frame, n_mut, n_wt = run_test(ids)
        if frame is None:
            log.append(f"  {subtype:<15}: skipped (mutant n={n_mut}, wild-type n={n_wt} - not enough lines)")
            continue
        log.append(f"  {subtype:<15}: mutant n={n_mut}, wild-type n={n_wt}")

        selected = frame[
            (frame["p_value"] < p_threshold)
            & (frame["difference"] <= min_effect_difference)
            & (frame["mutant_mean"] <= max_mutant_mean_effect)
        ]
        # "Masked" = significant inside this subtype, clearly not significant in the pooled test.
        masked = selected.join(pooled[["p_value", "difference"]], rsuffix="_pooled", how="inner")
        masked = masked[
            (masked["p_value_pooled"] >= p_threshold) | (masked["difference_pooled"] > min_effect_difference)
        ]
        masked = masked.join(pan_essential.rename("pct_all_dependent") * 100)
        masked = masked[masked["pct_all_dependent"] < 80]
        masked["subtype"] = subtype
        findings.append(masked.sort_values("difference"))

    if not findings or all(len(f) == 0 for f in findings):
        log.append("")
        log.append("No subtype-restricted candidate survived the filters.")
        return "\n".join(log)

    combined = pd.concat(findings).sort_values("difference")
    log.append("")
    log.append(f"STEP 3 | Candidates visible ONLY inside a subtype (top {min(top_n, len(combined))})")
    log.append(f"{'gene':<12}{'subtype':<16}{'MUT':>8}{'WT':>8}{'diff':>8}{'p':>10}{'p(pooled)':>11}{'%all dep':>10}")
    log.append("-" * 83)
    for gene, row in combined.head(top_n).iterrows():
        log.append(
            f"{gene:<12}{row['subtype']:<16}{row['mutant_mean']:>8.3f}{row['wildtype_mean']:>8.3f}"
            f"{row['difference']:>8.3f}{row['p_value']:>10.2e}{row['p_value_pooled']:>11.2e}"
            f"{row['pct_all_dependent']:>9.0f}%"
        )

    genes = list(dict.fromkeys(combined.head(top_n).index.tolist()))
    log.append("")
    log.append(f"SUBTYPE_MASKED_CANDIDATES: {', '.join(genes)}")
    log.append("")
    log.append("HOW TO READ THIS (important)")
    log.append("  These are NOT stronger candidates than the pooled hits - they are weaker in the sense that")
    log.append("  they rest on fewer cell lines, and subgroup testing inflates false positives. What makes")
    log.append("  them interesting is that the cell-line panel cannot in principle settle them: the relevant")
    log.append("  subtype is represented by a handful of lines. A subtype-matched organoid panel is the")
    log.append("  cheapest way to find out whether the signal is real.")
    log.append("  Next: assess_organoid_transferability(candidate_genes=[...]) to check whether organoid")
    log.append("  medium would mask the dependency before designing the experiment.")
    log.append("")
    log.append("PROVENANCE")
    for item in bundle["provenance"]:
        log.append(f"  - {item}")
    log.append(f"  - {expression_bundle['provenance']}")
    log.append(f"  - Mutation calls: {mutation_source}")
    log.append("  - No multiple-testing correction is applied to the subgroup tests (q-values would be")
    log.append("    uninformative at these group sizes); treat every row as a hypothesis.")
    return "\n".join(log)


def assess_organoid_transferability(
    candidate_genes: list[str] | str,
    target_mutation: str = "KRAS",
    cancer_type: str = "Pancreatic Cancer",
    normal_tissue: str = "Pancreas",
    data_lake_path: str | None = None,
) -> str:
    """Predict whether a cell-line synthetic lethal hit will still be detectable in an organoid.

    Organoid culture is not a more realistic version of the same experiment - it changes which genes
    are limiting. The medium supplies WNT3A/RSPO1, EGF, FGF10, Noggin, a TGF-beta receptor inhibitor
    and a ROCK inhibitor, so dependencies upstream of those factors vanish; and the cells grow
    anchored in laminin-rich matrix rather than on plastic, so adhesion and mechanotransduction genes
    matter differently. This tool checks each candidate against those changes and says what to expect.

    Four checks per candidate:
    * Niche coupling - is the gene upstream of a factor the medium supplies (dependency will be
      masked) or downstream of it (preserved)?
    * Anchorage sensitivity - measured directly as the gene-effect difference between adherent and
      suspension cell lines across the whole DepMap panel, a proxy for culture-geometry dependence.
    * Therapeutic window - expression in the matched normal tissue (GTEx), which is what a matched
      normal organoid arm would test.
    * Orthogonal genetic interaction - human genetic interactions reported in BioGRID, independent
      of DepMap.

    Parameters
    ----------
    candidate_genes : list[str] | str
        Candidate genes; a list, or a comma/whitespace separated string.
    target_mutation : str, optional
        The driver gene the candidates were called against (default: "KRAS").
    cancer_type : str, optional
        Cancer context, used only for reporting (default: "Pancreatic Cancer").
    normal_tissue : str, optional
        GTEx tissue used as the normal counterpart (default: "Pancreas").
    data_lake_path : str, optional
        Directory holding the DepMap and GTEx files.

    Returns
    -------
    str
        A research log with, per candidate, the niche-coupling verdict, the measured anchorage
        sensitivity, normal-tissue expression, orthogonal genetic-interaction support, an
        ORGANOID-ENHANCED / ORGANOID-MASKED / TRANSFERABLE / UNCERTAIN verdict, and the medium
        modification needed to keep the experiment interpretable.

    """
    from scipy import stats

    genes = _parse_gene_list(candidate_genes)
    if not genes:
        return "FAILURE: no candidate genes supplied."

    try:
        bundle = _load_depmap(data_lake_path)
    except FileNotFoundError as e:
        return f"FAILURE: {e}"

    gene_effect = bundle["gene_effect"]
    columns = bundle["gene_columns"]
    models = bundle["model"].set_index("ModelID")
    growth = models.reindex(gene_effect.index)["GrowthPattern"]
    adherent = gene_effect.index[growth == "Adherent"]
    suspension = gene_effect.index[growth == "Suspension"]

    gtex = _load_gtex(data_lake_path)
    interactions = _load_human_genetic_interactions(data_lake_path)

    log = [
        "=" * 78,
        f"ORGANOID TRANSFERABILITY - {target_mutation.upper()} SL candidates in {cancer_type}",
        f"Candidates: {', '.join(genes)}",
        f"Run at {datetime.now(tz=_UTC).strftime('%Y-%m-%d %H:%M UTC')}",
        "=" * 78,
        "",
        f"Anchorage reference panel: {len(adherent)} adherent vs {len(suspension)} suspension lines.",
        "Medium model: standard PDAC organoid feeding medium (WNT3A/RSPO1, EGF, FGF10, Noggin,",
        "A83-01, Y-27632) in laminin-rich matrix.",
    ]

    if gtex is None:
        log.append("GTEx table unavailable - the normal-tissue window check is skipped.")
    if interactions is None:
        log.append("Genetic-interaction table unavailable - the orthogonal check is skipped.")

    verdicts = {}
    for gene in genes:
        log.append("")
        log.append(f"### {gene}")
        if gene not in columns:
            log.append("  Not in the CRISPR library; skipped.")
            continue

        # --- niche coupling ---
        masked_by, preserved_by = [], []
        for factor, spec in NICHE_FACTORS.items():
            if gene in spec["masked_upstream"]:
                masked_by.append((factor, spec["note"]))
            if gene in spec["preserved_downstream"]:
                preserved_by.append(factor)
        if masked_by:
            for factor, note in masked_by:
                log.append(f"  Niche coupling: MASKED by '{factor}' - {note}")
        elif preserved_by:
            log.append(f"  Niche coupling: downstream of {', '.join(preserved_by)}; not rescued by the medium")
        else:
            log.append("  Niche coupling: no known overlap with a supplied niche factor")

        # --- anchorage sensitivity (measured) ---
        a = gene_effect.loc[adherent, columns[gene]].dropna()
        s = gene_effect.loc[suspension, columns[gene]].dropna()
        anchorage_difference, anchorage_p = float("nan"), 1.0
        if len(a) > 10 and len(s) > 10:
            anchorage_difference = a.mean() - s.mean()
            anchorage_p = stats.ttest_ind(a, s, equal_var=False)[1]
            if anchorage_p >= 0.01:
                direction = "no difference between anchored and suspension growth"
            elif anchorage_difference < 0:
                direction = "stronger when anchored"
            else:
                direction = "weaker when anchored"
            log.append(
                f"  Anchorage sensitivity: adherent {a.mean():+.3f} vs suspension {s.mean():+.3f} "
                f"(diff {anchorage_difference:+.3f}, p={anchorage_p:.1e}) - {direction}"
            )

        # --- therapeutic window via normal tissue ---
        normal_tpm = None
        if gtex is not None:
            row = gtex[(gtex.Tissue == normal_tissue) & (gtex.Gene == gene)]
            if len(row):
                normal_tpm = float(row.Expression.iloc[0])
                log.append(f"  Normal {normal_tissue} expression (GTEx): {normal_tpm:.1f} TPM")

        # --- orthogonal genetic interaction ---
        # Coverage is narrow (a few published combinatorial screens), so "no interaction found" must
        # be reported as "not screened", never as evidence that no interaction exists.
        if interactions is not None:
            driver = target_mutation.upper()
            covered = set(interactions.symbol_a) | set(interactions.symbol_b)
            if gene not in covered or driver not in covered:
                missing = [g for g in (gene, driver) if g not in covered]
                log.append(
                    f"  Orthogonal genetic interactions: {', '.join(missing)} not covered by the available "
                    f"human screens ({len(covered)} genes total) - absence of evidence, not evidence of absence"
                )
            else:
                pairs = interactions[
                    ((interactions.symbol_a == gene) & (interactions.symbol_b == driver))
                    | ((interactions.symbol_a == driver) & (interactions.symbol_b == gene))
                ]
                partners = len(interactions[(interactions.symbol_a == gene) | (interactions.symbol_b == gene)])
                log.append(
                    f"  Orthogonal genetic interactions (BioGRID human): {len(pairs)} with {driver}, "
                    f"{partners} partners in total"
                )

        # --- verdict ---
        if masked_by:
            factor = masked_by[0][0]
            verdict = (
                f"ORGANOID-MASKED - standard medium supplies {factor}, which rescues this dependency. "
                f"Drop {factor} from the medium (or use a reduced-factor formulation) or the organoid "
                "will return a false negative."
            )
        elif anchorage_p < 1e-3 and anchorage_difference < -0.15:
            verdict = (
                "ORGANOID-ENHANCED - the dependency is measurably stronger in anchored culture, so a "
                "3D matrix-embedded organoid should show a larger effect than the 2D screen did."
            )
        elif anchorage_p < 1e-3 and anchorage_difference > 0.15:
            verdict = (
                "ORGANOID-WEAKENED - the dependency is weaker in anchored culture; expect a smaller "
                "effect size in an organoid and power the experiment accordingly."
            )
        elif preserved_by:
            verdict = "TRANSFERABLE - downstream of the supplied niche factors and not anchorage-sensitive."
        else:
            verdict = (
                "UNCERTAIN - no medium conflict, and any anchorage effect is too small to change the "
                "readout; the organoid result should broadly track the cell-line result, so treat the "
                "organoid as a subtype and patient-diversity test rather than a mechanism test."
            )
        if normal_tpm is not None and normal_tpm >= 20:
            verdict += (
                f" WINDOW RISK: {normal_tpm:.0f} TPM in normal {normal_tissue}; run the matched normal "
                "organoid arm before any efficacy claim."
            )
        verdicts[gene] = verdict
        log.append(f"  VERDICT: {verdict}")

    log.append("")
    log.append("SUMMARY")
    for gene, verdict in verdicts.items():
        log.append(f"  {gene:<12}{verdict.split(' - ')[0]}")

    log.append("")
    log.append("PROVENANCE AND LIMITS")
    log.append(
        "  - No organoid dependency data is used or available: the local DepMap snapshot has 5 "
        "pancreatic organoid models (PANFR series) with NO CRISPR and NO expression data."
    )
    log.append(
        "  - Anchorage sensitivity is measured from adherent-vs-suspension cell lines, which is a "
        "proxy for culture geometry and is confounded with lineage (suspension lines are mostly "
        "haematopoietic)."
    )
    log.append(
        "  - Niche coupling is a curated knowledge map of standard organoid medium components "
        "(Boj et al., Cell 2015), not a measurement."
    )
    if gtex is not None:
        log.append(f"  - Normal tissue expression: GTEx median TPM, tissue '{normal_tissue}'.")
    if interactions is not None:
        covered_genes = len(set(interactions.symbol_a) | set(interactions.symbol_b))
        log.append(
            f"  - Genetic interactions: BioGRID human (taxid 9606), {len(interactions)} pairs covering only "
            f"{covered_genes} genes from a handful of combinatorial screens - most genes are simply untested."
        )
    return "\n".join(log)


def design_organoid_sl_experiment(
    candidate_gene: str,
    target_mutation: str = "KRAS",
    cancer_type: str = "Pancreatic Cancer",
    subtype: str = "both",
    data_lake_path: str | None = None,
) -> str:
    """Write a patient-derived organoid protocol that can confirm or kill a synthetic lethal candidate.

    Produces a concrete, criticism-ready plan: which organoid models to use and why, how to modify
    the medium given the candidate's niche coupling, how to deliver CRISPR into organoids, which 3D
    readouts to take, the matched normal-organoid arm that establishes the therapeutic window, and
    the pre-specified criteria that would falsify the hypothesis rather than confirm it.

    Parameters
    ----------
    candidate_gene : str
        The candidate synthetic lethal partner to test, e.g. "TEAD1".
    target_mutation : str, optional
        The driver gene defining the genotype contrast (default: "KRAS").
    cancer_type : str, optional
        Cancer context (default: "Pancreatic Cancer").
    subtype : str, optional
        Which transcriptional subtype the organoid panel should cover: "classical", "basal" or
        "both" (default: "both").
    data_lake_path : str, optional
        Directory holding the DepMap files, used to list organoid models present in DepMap.

    Returns
    -------
    str
        A step-by-step organoid experiment protocol with model selection, medium formulation,
        perturbation strategy, readouts, controls, powering and falsification criteria.

    """
    gene = candidate_gene.strip().upper()
    driver = target_mutation.strip().upper()

    masked_by = [(factor, spec["note"]) for factor, spec in NICHE_FACTORS.items() if gene in spec["masked_upstream"]]
    preserved_by = [factor for factor, spec in NICHE_FACTORS.items() if gene in spec["preserved_downstream"]]

    organoid_models = []
    try:
        bundle = _load_depmap(data_lake_path)
        organoids = bundle["model"]
        organoids = organoids[(organoids.get("ModelType") == "Organoid")]
        matched, _ = _select_cancer_models(organoids, cancer_type)
        organoid_models = matched["StrippedCellLineName"].dropna().tolist()
    except (FileNotFoundError, KeyError):
        pass

    log = [
        "=" * 78,
        f"ORGANOID VALIDATION PROTOCOL - {driver}-mutant {cancer_type} vs {gene} loss",
        f"Generated {datetime.now(tz=_UTC).strftime('%Y-%m-%d %H:%M UTC')}",
        "=" * 78,
        "",
        f"HYPOTHESIS: {gene} loss is selectively lethal in {driver}-mutant {cancer_type} organoids,",
        f"and this selectivity is not reproduced in {driver}-wild-type or normal organoids.",
        "",
        "1 | MODEL PANEL",
        f"   - {driver}-mutant tumour organoids: at least 6 independent patient lines.",
    ]
    if subtype == "both":
        log.append("   - Subtype coverage: >=3 classical and >=3 basal-like lines, typed by GATA6/KRT81")
        log.append("     expression before the experiment. Subtype is the variable a cell-line panel could")
        log.append("     not control, so losing it here wastes the model system.")
    else:
        log.append(f"   - Subtype coverage: {subtype}-restricted panel, typed by GATA6/KRT81 before use.")
    log.append(f"   - {driver}-wild-type tumour organoids: >=3 lines (rare; source deliberately).")
    log.append("   - Matched NORMAL organoids: >=3 lines from adjacent non-tumour tissue of the same")
    log.append("     patients where possible. This is the arm that cell lines cannot provide.")
    if organoid_models:
        log.append(f"   - Organoid models present in DepMap for this context: {', '.join(organoid_models)}")
        log.append("     (annotation only - DepMap has no CRISPR or expression data for them, so they")
        log.append("      must be screened in-house or sourced from a biobank such as HCMI.)")

    log.append("")
    log.append("2 | NICHE SUBTYPE TYPING (do this BEFORE choosing a medium)")
    log.append("   PDAC organoids are not uniform in what the niche must supply. Seino et al.")
    log.append("   (Cell Stem Cell 2018, PMID 29337182) typed 39 patient-derived PDAC organoids into")
    log.append("   three functional subtypes by Wnt/R-spondin dependency:")
    log.append("     (a) Wnt-non-producing  - requires Wnt from cancer-associated fibroblasts")
    log.append("     (b) Wnt-producing      - secretes its own Wnt ligands")
    log.append("     (c) R-spondin-independent - grows without Wnt and R-spondin")
    log.append("   These subtypes track GATA6-dependent expression subtypes, i.e. the same")
    log.append("   classical/basal axis used in discover_subtype_masked_sl_candidates.")
    log.append("   Procedure: grow each line in complete / Wnt-withdrawn / Wnt+RSPO-withdrawn medium")
    log.append("   for 2 passages and record which conditions sustain growth.")
    log.append("   Why it matters here: feeding every line the same WNT3A/RSPO1-replete medium erases")
    log.append("   the very heterogeneity that distinguishes patients, and a knockout phenotype seen")
    log.append("   only in Wnt-non-producing lines will look irreproducible if the subtype is ignored.")
    log.append("")
    log.append("3 | MEDIUM")
    log.append("   Base: human complete feeding medium (WNT3A/RSPO1, Noggin, EGF, FGF10, A83-01,")
    log.append("   nicotinamide, gastrin, N-acetylcysteine, B27) in laminin-rich matrix domes.")
    log.append("   Match the base to the niche subtype determined in step 2; report the subtype")
    log.append("   composition of the panel alongside any result.")
    if masked_by:
        for factor, note in masked_by:
            log.append(f"   !! CRITICAL MODIFICATION: {gene} is upstream of '{factor}' ({note}).")
            log.append("      Standard medium will rescue the knockout and produce a FALSE NEGATIVE.")
            log.append(f"      Run paired arms: complete medium vs medium WITHOUT {factor}.")
            log.append("      A phenotype that appears only in the withdrawal arm is still a real")
            log.append("      dependency, but it is conditional - report it that way.")
    elif preserved_by:
        log.append(f"   {gene} acts downstream of {', '.join(preserved_by)}; standard medium is appropriate.")
    else:
        log.append("   No known conflict between this gene and the supplied niche factors.")
    log.append("   Y-27632 only during seeding (48 h); it masks ROCK-dependent phenotypes if kept on.")

    log.append("")
    log.append("4 | PERTURBATION")
    log.append(f"   - Doxycycline-inducible Cas9 organoid lines; two independent {gene} sgRNAs plus")
    log.append("     non-targeting and an essential-gene (e.g. PLK1) positive control.")
    log.append("   - Transduce dissociated single cells, select, then re-embed; confirm editing by")
    log.append("     amplicon sequencing and protein loss by western blot in 3D lysates.")
    log.append("   - Induce AFTER organoids are established, so the readout is maintenance rather than")
    log.append("     formation - otherwise seeding efficiency confounds everything.")

    log.append("")
    log.append("5 | READOUTS (3D-appropriate)")
    log.append("   - Organoid formation efficiency and size distribution (bright-field, automated).")
    log.append("   - 3D CellTiter-Glo at day 7 and day 14.")
    log.append("   - Cleaved caspase-3 by whole-mount immunofluorescence to separate death from arrest.")
    log.append("   - Passage capacity over 3 serial passages: the phenotype that matters clinically is")
    log.append("     loss of regrowth, not a day-7 viability dip.")

    log.append("")
    log.append("6 | POWERING")
    log.append("   - Unit of replication is the PATIENT LINE, not the well. n=6 mutant vs n=3 wild-type")
    log.append("     lines detects a 2-fold selective effect at ~80% power only if between-line variance")
    log.append("     is modest; run >=3 technical replicates per line and analyse with a mixed model")
    log.append("     (line as random effect, genotype x sgRNA as fixed effects).")

    log.append("")
    log.append("7 | FALSIFICATION CRITERIA (pre-specified)")
    log.append("   The hypothesis is REJECTED if any of these holds:")
    log.append(f"   - Normal organoids show comparable {gene}-loss sensitivity -> no therapeutic window.")
    log.append("   - Only one of the two sgRNAs produces the phenotype -> off-target.")
    log.append(f"   - Re-expression of sgRNA-resistant {gene} fails to rescue -> off-target.")
    log.append(f"   - The effect is present in {driver}-wild-type organoids at similar magnitude ->")
    log.append("     a general dependency, not synthetic lethality.")
    log.append("   - The effect appears only in one subtype -> report as subtype-restricted, and do not")
    log.append("     generalise to the disease.")
    if masked_by:
        log.append("   - The effect appears ONLY in the factor-withdrawal arm -> conditional dependency;")
        log.append("     it does not support targeting the gene in a niche-replete tumour.")

    log.append("")
    log.append("8 | WHAT THIS EXPERIMENT CANNOT SETTLE")
    log.append("   - Stromal and immune contributions: tumour organoid monoculture has neither.")
    log.append("     A negative result here does not exclude a dependency that requires the stroma.")
    log.append("   - Pharmacology: genetic loss is not the same as inhibition of the encoded protein.")
    log.append("   - Patient-level efficacy: organoid response correlates with, but does not establish,")
    log.append("     clinical benefit.")
    return "\n".join(log)
