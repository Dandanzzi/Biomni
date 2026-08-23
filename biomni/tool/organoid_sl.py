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
3. ``discover_serum_masked_dependencies`` - organoid feeding medium is serum-free, so cholesterol
   and fatty acids must be made de novo, while a 2D screen in 10% FBS hands the cell its sterols.
   Whole classes of dependency are therefore invisible to the cell-line panel by construction. This
   tool measures the effect on the cell lines DepMap does annotate as serum-free, lineage-matched.
4. ``discover_allele_resolved_sl_candidates`` - "KRAS-mutant" is not one genotype. Pooling G12D,
   G12V, G12R, G12C and Q61H into a single label averages away allele-restricted dependencies. This
   tool repeats the test one allele at a time and reports what the pooled label hid.
5. ``design_organoid_sl_experiment`` - turns the above into a concrete PDO protocol with a matched
   normal-organoid arm, which is the therapeutic-window control that no cell-line panel provides.

Blind spots 3 and 4 follow the design of the tumour-derived organoid biobank study (Nature, 2026;
256 patient-derived organoids, genome-wide CRISPR in 162 of them). That study reported 97
organoid-specific core fitness genes on top of the 654 shared with the cell-line DepMap, enriched
for steroid / cholesterol / isoprenoid biosynthesis, and showed that KRAS dependency in colorectal
organoids splits by allele (G12 alleles retain KRAS and EGFR/PTPN11 dependency; Q61H does not).
Both findings are properties of the *model system and the genotype label*, not of one cancer type,
so both can be probed on the local DepMap snapshot before committing to an organoid screen.

Honest limitation, stated up front and repeated in every tool's output: the local DepMap snapshot
contains 5 pancreatic organoid models (PANFR series) but **no CRISPR or expression data for them**
(0 of the 24 organoid models in the snapshot carry either). Nothing here measures an organoid
dependency. These tools predict where cell-line evidence is untrustworthy and design the organoid
experiment that would settle it.
"""

import os
import re
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
_SERUM_CONTRAST_CACHE: dict = {}

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
    "Serum-free base (B27/N2, no lipoprotein)": {
        # This entry runs the opposite way to the others: the medium WITHHOLDS something that a
        # standard 2D screen supplies. Uptake genes are the ones rendered pointless. `withholds`
        # tells the verdict text to say "add it back", not "leave it out".
        "withholds": True,
        "masked_upstream": ["LDLR", "PCSK9", "MYLIP", "SCARB1", "LRP1", "NPC1", "NPC2", "LIPA", "CD36"],
        "preserved_downstream": [],
        "note": (
            "organoid medium carries no serum lipoprotein, so cholesterol uptake genes are dispensable "
            "while the de novo sterol pathway becomes limiting - the inverse of a 10% FBS screen"
        ),
    },
}

# Organoid-specific core fitness programmes (Nature 2026 organoid biobank: 97 organoid-specific core
# fitness genes on top of the 654 shared with the cell-line DepMap, enriched for steroid / cholesterol
# / isoprenoid biosynthesis).
#
# Why this is measurable here rather than merely asserted: the enrichment has a mechanical cause that
# DepMap annotates. Organoid feeding medium is serum-free (B27/N2 base), so there is no lipoprotein
# source and sterols must be synthesised; a standard cell-line screen runs in 10% FBS and supplies
# them. DepMap flags 6 CRISPR-screened lines as SerumFreeMedia, which makes the same contrast
# testable on cell lines - see `discover_serum_masked_dependencies`.
#
# The programmes are split rather than lumped because they do not behave alike: the sensing/processing
# and post-squalene arms are the ones that respond to sterol starvation, while the upper mevalonate
# steps are shared with the non-sterol isoprenoid branch and need not follow.
SERUM_DEPENDENT_PROGRAMS = {
    "SREBP sensing and processing": ["SCAP", "MBTPS1", "MBTPS2", "SREBF1", "SREBF2", "INSIG1", "INSIG2"],
    "cholesterol synthesis (post-squalene)": [
        "FDFT1",
        "SQLE",
        "LSS",
        "CYP51A1",
        "TM7SF2",
        "MSMO1",
        "NSDHL",
        "HSD17B7",
        "EBP",
        "SC5D",
        "DHCR7",
        "DHCR24",
    ],
    "mevalonate and isoprenoid backbone": [
        "ACAT2",
        "HMGCS1",
        "HMGCR",
        "MVK",
        "PMVK",
        "MVD",
        "IDI1",
        "FDPS",
        "GGPS1",
    ],
    "de novo lipogenesis": ["ACLY", "ACACA", "FASN", "SCD", "ELOVL1", "ELOVL5", "ELOVL6", "ACSL3"],
    "lipoprotein uptake (inverse control)": [
        "LDLR",
        "PCSK9",
        "MYLIP",
        "SCARB1",
        "LRP1",
        "NPC1",
        "NPC2",
        "LIPA",
        "CD36",
    ],
}

# Read-outs for the allele-resolved contrast. The organoid biobank study found that colorectal
# organoids carrying KRAS G12 alleles kept a strong KRAS self-dependency together with EGFR and
# PTPN11 dependency, whereas KRAS Q61H organoids were refractory both to EGFR inhibition and to EGF
# withdrawal. Reporting these genes per allele is a direct check of whether the same split exists in
# another lineage.
RAS_PATHWAY_PROBES = [
    "KRAS",
    "EGFR",
    "PTPN11",
    "SOS1",
    "SOS2",
    "SHOC2",
    "RAF1",
    "BRAF",
    "MAP2K1",
    "MAPK1",
    "RASA1",
    "SPRED2",
    "NRAS",
    "HRAS",
]


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


_MISSENSE_RE = re.compile(r"^p?\.?([A-Z])(\d+)([A-Z*])$")


def _normalise_allele(variant: str) -> tuple[str, str]:
    """Turn a protein-change string into a (allele, codon) pair, e.g. "p.G12D" -> ("G12D", "G12").

    A line carrying several distinct protein changes in the same gene is labelled MULTI rather than
    being silently assigned to one of them; a non-missense call becomes OTHER. Both are reported and
    excluded from allele-level testing instead of being folded into a neighbouring group.
    """
    tokens = [t.strip() for t in str(variant).split(",") if t.strip()]
    hits = {}
    for token in tokens:
        match = _MISSENSE_RE.match(token)
        if match:
            hits[f"{match.group(1)}{match.group(2)}{match.group(3)}"] = f"{match.group(1)}{match.group(2)}"
    if len(hits) == 1:
        allele, codon = next(iter(hits.items()))
        return allele, codon
    if len(hits) > 1:
        return "MULTI(" + "/".join(sorted(hits)) + ")", "MULTI"
    return "OTHER", "OTHER"


def _permutation_p(values, subset_genes, n_permutations: int, direction: str = "less", seed: int = 0):
    """Empirical p-value for the mean of ``subset_genes`` against random gene sets of equal size.

    Used instead of a parametric test because the per-gene statistics being aggregated are z-scores
    from only a handful of lines: their null distribution is not something to assume. A fixed seed
    keeps the whole pipeline reproducible.
    """
    import numpy as np

    present = [g for g in subset_genes if g in values.index]
    if len(present) < 3:
        return float("nan"), float("nan"), present
    observed = float(values[present].mean())
    pool = values.to_numpy(dtype=float)
    rng = np.random.default_rng(seed)
    draws = np.empty(n_permutations, dtype=float)
    for i in range(n_permutations):
        draws[i] = rng.choice(pool, size=len(present), replace=False).mean()
    if direction == "less":
        pvalue = float((draws <= observed).mean())
    else:
        pvalue = float((draws >= observed).mean())
    # An empirical p of exactly 0 only means "smaller than 1/n_permutations".
    pvalue = max(pvalue, 1.0 / n_permutations)
    return observed, pvalue, present


def _serum_free_contrast(data_lake_path: str | None = None, min_peers: int = 8):
    """Lineage-matched gene-effect z-scores for the DepMap lines that are cultured without serum.

    Organoid feeding medium is serum-free, so cholesterol and fatty acids have to be made rather
    than taken up. DepMap annotates a small number of CRISPR-screened lines as ``SerumFreeMedia``,
    which makes that contrast measurable on cell lines - but comparing them straight against the
    rest of the panel would confound medium with lineage and growth pattern, since the serum-free
    lines are not a random sample of the panel. Each serum-free line is therefore scored against
    serum-cultured lines of its OWN lineage (and growth pattern where enough exist), and the
    per-gene z-scores are averaged across lines.

    Returns None if the annotation is absent or no serum-free line has enough matched peers.
    """
    import numpy as np
    import pandas as pd

    resolved = _resolve_data_lake(data_lake_path)
    cache_key = (resolved, min_peers)
    if cache_key in _SERUM_CONTRAST_CACHE:
        return _SERUM_CONTRAST_CACHE[cache_key]

    bundle = _load_depmap(data_lake_path)
    gene_effect = bundle["gene_effect"]
    models = bundle["model"].set_index("ModelID").reindex(gene_effect.index)
    if "SerumFreeMedia" not in models.columns:
        return None

    serum_free_ids = models.index[models["SerumFreeMedia"] == True].tolist()  # noqa: E712
    serum_mask = models["SerumFreeMedia"] == False  # noqa: E712
    if not serum_free_ids or not serum_mask.any():
        return None

    lineage = models["OncotreeLineage"]
    growth = models["GrowthPattern"]

    z_by_model, rows = {}, []
    for model_id in serum_free_ids:
        peer_mask = serum_mask & (lineage == lineage.loc[model_id]) & (growth == growth.loc[model_id])
        matched_on = "lineage + growth pattern"
        if int(peer_mask.sum()) < min_peers:
            peer_mask = serum_mask & (lineage == lineage.loc[model_id])
            matched_on = "lineage only"
        n_peers = int(peer_mask.sum())
        if n_peers < min_peers:
            rows.append(
                {
                    "model": models.loc[model_id, "StrippedCellLineName"],
                    "lineage": lineage.loc[model_id],
                    "growth": growth.loc[model_id],
                    "n_peers": n_peers,
                    "matched_on": "SKIPPED - too few serum-cultured peers",
                }
            )
            continue
        block = gene_effect.loc[peer_mask[peer_mask].index]
        z_by_model[model_id] = (gene_effect.loc[model_id] - block.mean()) / block.std().replace(0, np.nan)
        rows.append(
            {
                "model": models.loc[model_id, "StrippedCellLineName"],
                "lineage": lineage.loc[model_id],
                "growth": growth.loc[model_id],
                "n_peers": n_peers,
                "matched_on": matched_on,
            }
        )

    if not z_by_model:
        return None

    z_frame = pd.DataFrame(z_by_model).T
    mean_z = z_frame.mean()
    mean_z.index = [c.split(" (")[0].upper() for c in gene_effect.columns]
    n_scored = z_frame.notna().sum()
    n_scored.index = mean_z.index
    result = {
        "z": mean_z.dropna(),
        "model_ids": list(z_by_model),
        "n_models_scored": n_scored,
        "models": pd.DataFrame(rows),
        "n_serum_free": len(z_by_model),
        "n_serum": int(serum_mask.sum()),
        "provenance": bundle["provenance"],
    }
    _SERUM_CONTRAST_CACHE[cache_key] = result
    return result


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
                "gene": [c.split(" (")[0].upper() for c in usable],
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
    pan_essential.index = [c.split(" (")[0].upper() for c in gene_effect.columns]

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


def discover_serum_masked_dependencies(
    cancer_type: str = "Pancreatic Cancer",
    target_mutation: str | None = "KRAS",
    data_lake_path: str | None = None,
    mutation_csv_path: str | None = None,
    n_permutations: int = 5000,
    min_models_scored: int = 4,
    z_threshold: float = -1.0,
    max_serum_free_effect: float = -0.5,
    min_cohort_expression: float = 1.0,
    top_n: int = 20,
    max_pct_all_dependent: float = 80.0,
    output_csv_path: str | None = None,
) -> str:
    """Find dependencies that a serum-grown cell-line screen cannot see but a serum-free organoid can.

    The tumour-derived organoid biobank study (Nature, 2026) screened 162 patient-derived organoids
    genome-wide and found 97 core fitness genes that are organoid-specific - absent from the
    cell-line dependency map - enriched for steroid, cholesterol and isoprenoid biosynthesis. That
    enrichment is not a mystery of three-dimensionality: organoid feeding medium is serum-free, so
    the cells must synthesise their own sterols and fatty acids, whereas a standard screen in 10%
    foetal bovine serum supplies them. A dependency that only exists when the nutrient is absent is
    invisible to the cell-line panel by construction, not by chance.

    DepMap annotates a small number of CRISPR-screened lines as serum-free, so the contrast can be
    measured rather than assumed. Each serum-free line is z-scored against serum-cultured lines of
    its own lineage and growth pattern, so the comparison is not simply reading out lineage. The
    tool then reports (a) whether the lipid programmes are shifted as a block, tested by permutation
    against random gene sets, and (b) which individual genes look non-essential in the cancer panel
    yet are strongly serum-masked - the candidates a 2D screen of this tumour type would have thrown
    away.

    Parameters
    ----------
    cancer_type : str, optional
        Cancer context whose 2D dependency is used as the "would we have found it?" reference
        (default: "Pancreatic Cancer").
    target_mutation : str, optional
        Driver gene; when given, the mutant-versus-wild-type difference in that cohort is reported
        alongside each candidate. Pass None to skip genotype stratification (default: "KRAS").
    data_lake_path : str, optional
        Directory holding the DepMap files (default: "./data/biomni_data/data_lake").
    mutation_csv_path : str, optional
        Optional mutation table overriding the default mutation source.
    n_permutations : int, optional
        Permutations used for the programme-level enrichment test (default: 5000).
    min_models_scored : int, optional
        A gene must be measured in at least this many serum-free lines to be ranked (default: 4).
    z_threshold : float, optional
        Lineage-matched z-score at or below which a gene counts as serum-masked (default: -1.0).
    max_serum_free_effect : float, optional
        The gene must also be a real dependency in the serum-free lines themselves: its mean gene
        effect there must be at or below this value (default: -0.5). Without this, the ranking fills
        up with genes whose effect is ~0 everywhere and whose z-score is pure noise.
    min_cohort_expression : float, optional
        Minimum mean expression in the cancer cohort, log2(TPM+1) (default: 1.0). Removes genes that
        are not expressed in this tumour type, where a CRISPR score is usually a copy-number or
        multi-mapping artefact rather than a dependency.
    top_n : int, optional
        Number of serum-masked candidates to report (default: 20).
    max_pct_all_dependent : float, optional
        Genes essential in more than this percentage of all screened lines are dropped as
        pan-essential (default: 80.0).
    output_csv_path : str, optional
        Write the full ranked candidate table to this path as well.

    Returns
    -------
    str
        A research log with the serum-free panel and how each line was matched, a permutation test
        per lipid programme including the lipoprotein-uptake inverse control, a ranked table of
        serum-masked genes annotated with their dependency in the cancer cohort, and a
        MISSED-BY-2D flag marking the ones the cell-line panel would have discarded.

    """
    import pandas as pd

    log = [
        "=" * 78,
        f"SERUM-MASKED DEPENDENCY DISCOVERY - organoid axis, reference cohort {cancer_type}",
        f"Run at {datetime.now(tz=_UTC).strftime('%Y-%m-%d %H:%M UTC')}",
        "=" * 78,
        "",
        "Premise: organoid medium is serum-free, so sterols and fatty acids must be made de novo.",
        "A screen in 10% FBS hands the cell those molecules, so the entire biosynthetic arm can look",
        "dispensable. This is the mechanism behind the organoid-specific core fitness genes reported",
        "by the Nature 2026 organoid biobank (97 genes, enriched for steroid/cholesterol/isoprenoid",
        "biosynthesis), and DepMap annotates enough serum-free lines to measure it directly.",
        "",
    ]

    try:
        bundle = _load_depmap(data_lake_path)
    except FileNotFoundError as e:
        return f"FAILURE: {e}"

    contrast = _serum_free_contrast(data_lake_path)
    if contrast is None:
        return (
            "FAILURE: the DepMap model table has no usable SerumFreeMedia annotation, or no serum-free "
            "line has enough lineage-matched peers. Without it the serum axis cannot be measured."
        )

    log.append("STEP 1 | The serum-free panel, and how each line was matched")
    log.append(f"  {'line':<12}{'lineage':<20}{'growth':<12}{'peers':>7}  matched on")
    log.append("  " + "-" * 74)
    for _, row in contrast["models"].iterrows():
        log.append(
            f"  {str(row['model']):<12}{str(row['lineage']):<20}{str(row['growth']):<12}"
            f"{int(row['n_peers']):>7}  {row['matched_on']}"
        )
    log.append(
        f"  Usable serum-free lines: {contrast['n_serum_free']} against {contrast['n_serum']} serum-cultured lines."
    )
    log.append("  READ THIS BEFORE THE NUMBERS: this panel is tiny, none of these lines is from the")
    log.append(f"  {cancer_type.lower()} lineage, and they were not cultured serum-free as an experiment - it")
    log.append("  is simply how each line is maintained. The lineage matching removes the most obvious")
    log.append("  confounder but cannot remove the others. Everything below is a hypothesis about which")
    log.append("  dependencies the organoid medium will unmask, not a measurement of an organoid.")

    gene_effect = bundle["gene_effect"]
    z = contrast["z"]
    scored = contrast["n_models_scored"]
    z = z[scored.reindex(z.index).fillna(0) >= min_models_scored]

    log.append("")
    log.append(f"STEP 2 | Are the lipid programmes shifted as a block? ({n_permutations} permutations)")
    log.append(f"  Genome-wide mean lineage-matched z = {z.mean():+.3f} (this is the null the programmes are")
    log.append("  compared against; a negative z means a stronger dependency without serum)")
    log.append("")
    log.append(f"  {'programme':<42}{'n':>4}{'mean z':>9}{'perm p':>9}  verdict")
    log.append("  " + "-" * 74)
    program_rows = []
    for name, genes in SERUM_DEPENDENT_PROGRAMS.items():
        inverse = "inverse control" in name
        observed, pvalue, present = _permutation_p(z, genes, n_permutations, direction="greater" if inverse else "less")
        if present is None or len(present) < 3:
            log.append(f"  {name:<42}{'-':>4}{'-':>9}{'-':>9}  too few genes measured")
            continue
        if inverse:
            verdict = "dispensable without serum, as predicted" if pvalue < 0.05 else "no shift detected"
        else:
            verdict = "ENRICHED - unmasked without serum" if pvalue < 0.05 else "no shift detected"
        log.append(f"  {name:<42}{len(present):>4}{observed:>+9.3f}{pvalue:>9.4f}  {verdict}")
        program_rows.append({"programme": name, "n_genes": len(present), "mean_z": observed, "perm_p": pvalue})
    log.append("")
    log.append("  The uptake row is the internal control and is tested in the opposite direction: with no")
    log.append("  lipoprotein in the medium there is nothing to import, so those genes should become LESS")
    log.append("  important. If the biosynthetic programmes and the uptake programme moved the same way,")
    log.append("  the signal would be a culture artefact rather than sterol starvation.")

    # --- reference dependency in the cancer cohort ------------------------------------------------
    cohort, match_note = _select_cancer_models(bundle["model"], cancer_type)
    cohort = cohort[cohort["ModelID"].isin(gene_effect.index)]
    if len(cohort) < 5:
        return "\n".join(log + ["", f"FAILURE: only {len(cohort)} screened lines matched '{cancer_type}'."])

    cohort_effect = gene_effect.loc[cohort["ModelID"]]
    cohort_mean = cohort_effect.mean()
    cohort_pct = (cohort_effect < DEPLETION_THRESHOLD).sum() / cohort_effect.notna().sum() * 100
    pan_pct = (gene_effect < DEPLETION_THRESHOLD).sum() / gene_effect.notna().sum() * 100
    for series in (cohort_mean, cohort_pct, pan_pct):
        series.index = [c.split(" (")[0].upper() for c in gene_effect.columns]

    mutant_difference = None
    mutation_source = None
    if target_mutation:
        try:
            annotated, mutation_source = _annotate_mutation_status(
                cohort, target_mutation, bundle["data_lake_path"], mutation_csv_path
            )
            mutant_ids = annotated.loc[annotated["MutationStatus"] == "MUT", "ModelID"].tolist()
            wildtype_ids = annotated.loc[annotated["MutationStatus"] == "WT", "ModelID"].tolist()
            if len(mutant_ids) >= 3 and len(wildtype_ids) >= 2:
                mutant_difference = gene_effect.loc[mutant_ids].mean() - gene_effect.loc[wildtype_ids].mean()
                mutant_difference.index = [c.split(" (")[0].upper() for c in gene_effect.columns]
        except (RuntimeError, ValueError) as e:
            log.append("")
            log.append(f"  (genotype stratification skipped: {e})")

    serum_free_effect = gene_effect.loc[contrast["model_ids"]].mean()
    serum_free_effect.index = [c.split(" (")[0].upper() for c in gene_effect.columns]

    cohort_expression = None
    try:
        expression_bundle = _load_expression(data_lake_path)
        expressed_ids = [m for m in cohort["ModelID"] if m in expression_bundle["expression"].index]
        if len(expressed_ids) >= 5:
            cohort_expression = expression_bundle["expression"].loc[expressed_ids].mean()
            cohort_expression.index = [c.split(" (")[0].upper() for c in expression_bundle["expression"].columns]
    except FileNotFoundError:
        cohort_expression = None

    program_membership = {}
    for name, genes in SERUM_DEPENDENT_PROGRAMS.items():
        for gene in genes:
            program_membership.setdefault(gene, name)

    table = pd.DataFrame(
        {
            "gene": z.index,
            "serum_free_z": z.to_numpy(),
            "n_serum_free_scored": scored.reindex(z.index).to_numpy(),
            "cohort_mean_effect": cohort_mean.reindex(z.index).to_numpy(),
            "pct_cohort_dependent": cohort_pct.reindex(z.index).to_numpy(),
            "pct_all_dependent": pan_pct.reindex(z.index).to_numpy(),
        }
    )
    table["lipid_programme"] = table["gene"].map(program_membership).fillna("-")
    if mutant_difference is not None:
        table["mutant_minus_wt"] = mutant_difference.reindex(z.index).to_numpy()
    table["serum_free_mean_effect"] = serum_free_effect.reindex(z.index).to_numpy()
    if cohort_expression is not None:
        table["cohort_expression"] = cohort_expression.reindex(z.index).to_numpy()

    before = len(table)
    table = table[table["pct_all_dependent"] < max_pct_all_dependent]
    dropped_pan = before - len(table)
    before = len(table)
    # A z-score says "more depleted than its lineage peers"; it does not say the gene is a
    # dependency at all. Without this filter the ranking is topped by genes whose gene effect is ~0
    # in every line, where the z-score is noise divided by a small standard deviation.
    table = table[table["serum_free_mean_effect"] <= max_serum_free_effect]
    dropped_weak = before - len(table)
    dropped_unexpressed = 0
    if cohort_expression is not None:
        before = len(table)
        table = table[table["cohort_expression"] >= min_cohort_expression]
        dropped_unexpressed = before - len(table)
    # MISSED-BY-2D: the cancer panel, grown in serum, would not have called this gene a dependency.
    table["missed_by_2d"] = (table["cohort_mean_effect"] > DEPLETION_THRESHOLD) & (table["serum_free_z"] <= z_threshold)
    table = table.sort_values("serum_free_z")

    candidates = table[table["serum_free_z"] <= z_threshold]
    log.append("")
    log.append(f"STEP 3 | Genes most strongly unmasked without serum (top {min(top_n, len(candidates))})")
    log.append(f"  {len(candidates)} genes reach z <= {z_threshold} after three filters:")
    log.append(f"    - pan-essential (>= {max_pct_all_dependent:.0f}% of all lines dependent): {dropped_pan} dropped")
    log.append(
        f"    - not a dependency in the serum-free lines either (effect > {max_serum_free_effect}): "
        f"{dropped_weak} dropped"
    )
    if cohort_expression is not None:
        log.append(
            f"    - not expressed in the cohort (< {min_cohort_expression} log2 TPM+1): {dropped_unexpressed} dropped"
        )
    else:
        log.append("    - expression filter skipped: no expression matrix available")
    log.append("")
    header = f"  {'gene':<12}{'z':>7}{'n':>3}{'sf eff':>8}{'cohort eff':>12}{'%coh dep':>10}{'%all dep':>10}"
    if mutant_difference is not None:
        header += f"{'MUT-WT':>9}"
    header += "  programme / flag"
    log.append(header)
    log.append("  " + "-" * (len(header) + 10))
    for _, row in candidates.head(top_n).iterrows():
        line = (
            f"  {row['gene']:<12}{row['serum_free_z']:>7.2f}{int(row['n_serum_free_scored']):>3}"
            f"{row['serum_free_mean_effect']:>8.2f}"
            f"{row['cohort_mean_effect']:>12.3f}{row['pct_cohort_dependent']:>9.0f}%{row['pct_all_dependent']:>9.0f}%"
        )
        if mutant_difference is not None:
            line += f"{row['mutant_minus_wt']:>9.3f}"
        tag = row["lipid_programme"] if row["lipid_programme"] != "-" else ""
        if row["missed_by_2d"]:
            tag = (tag + " | MISSED-BY-2D").strip(" |")
        log.append(line + f"  {tag}")

    missed = candidates[candidates["missed_by_2d"]]
    log.append("")
    log.append("STEP 4 | The candidates worth an organoid")
    log.append(
        f"  {len(missed)} of the {len(candidates)} serum-masked genes are NOT dependencies in the "
        f"{cancer_type.lower()} panel"
    )
    log.append(
        f"  (cohort mean gene effect above {DEPLETION_THRESHOLD}), so a 2D screen of this tumour type "
        "would have discarded them."
    )
    log.append("  Those are the ones where the organoid is not a confirmation step but the only experiment")
    log.append("  that can give an answer.")
    reported = missed.head(top_n)["gene"].tolist()
    log.append("")
    log.append(f"SERUM_MASKED_CANDIDATES: {', '.join(reported) if reported else '(none)'}")

    if output_csv_path:
        table.to_csv(output_csv_path, index=False)
        log.append("")
        log.append(f"Full ranked table written to {output_csv_path} ({len(table)} genes).")

    log.append("")
    log.append("HOW TO READ THIS (important)")
    log.append("  A serum-masked gene is a prediction about the MEDIUM, not about three-dimensionality.")
    log.append("  It says: if you run this screen in organoid medium, this gene should score, and the")
    log.append("  reason your cell-line screen missed it is that the serum was feeding the cells.")
    log.append("  The direct falsification is cheap and does not need an organoid: re-run the candidate")
    log.append("  in the same cell lines in lipid-depleted serum, and the dependency should appear.")
    log.append("  Do that before spending an organoid screen on it.")
    log.append("  Note also that the therapeutic-window argument is the reverse of the usual one here:")
    log.append("  normal tissue in a patient is not lipid-starved, so a dependency that only exists")
    log.append("  without lipoprotein may be an artefact of the culture rather than a drug target.")
    log.append("  assess_organoid_transferability() reports the normal-tissue expression for each gene.")
    log.append("")
    log.append("PROVENANCE AND LIMITS")
    for item in bundle["provenance"]:
        log.append(f"  - {item}")
    log.append(f"  - Cohort: {match_note}")
    if mutation_source:
        log.append(f"  - Mutation calls: {mutation_source}")
    log.append(
        f"  - Serum axis measured on {contrast['n_serum_free']} serum-free lines only. Any single "
        "outlying line can move a gene several z-units; treat per-gene values as noisy and the "
        "programme-level permutation test as the more trustworthy statistic."
    )
    log.append(
        "  - No organoid CRISPR data is used: 0 of the 24 organoid models in the local DepMap snapshot "
        "have CRISPR or expression data."
    )
    log.append(
        "  - Serum-free medium is a proxy for one property of organoid medium. It says nothing about "
        "WNT/RSPO1, EGF, FGF10, Noggin, A83-01 or matrix embedding - assess_organoid_transferability() "
        "covers those."
    )
    return "\n".join(log)


def discover_allele_resolved_sl_candidates(
    cancer_type: str = "Pancreatic Cancer",
    target_mutation: str = "KRAS",
    data_lake_path: str | None = None,
    mutation_csv_path: str | None = None,
    min_allele_lines: int = 4,
    p_threshold: float = 0.05,
    min_effect_difference: float = -0.2,
    max_mutant_mean_effect: float = -0.3,
    max_pct_all_dependent: float = 80.0,
    top_n: int = 15,
    output_csv_path: str | None = None,
) -> str:
    """Find synthetic lethal candidates that the label "<gene>-mutant" hides by pooling alleles.

    "KRAS-mutant" is a label, not a genotype. The organoid biobank study (Nature, 2026) showed the
    cost of treating it as one: colorectal organoids carrying KRAS G12 alleles retained a strong
    KRAS self-dependency together with EGFR and PTPN11 dependency, while KRAS Q61H organoids were
    refractory both to EGFR inhibition and to EGF withdrawal. A dependency restricted to one allele
    is diluted by every line carrying a different allele, so a pooled mutant-versus-wild-type test
    can only find it if the allele happens to dominate the panel.

    This tool splits the mutant lines by protein change and repeats the test one allele at a time.
    It runs two contrasts, because they fail in different ways:

    * allele versus wild type - the same question the pooled test asks, but restricted. Limited by
      how few wild-type lines exist for a driver that is nearly universal in this tumour type.
    * allele versus the other mutant alleles - the direct test of whether the pooled label is
      hiding something. It needs no wild-type lines at all, which is what makes it usable when the
      wild-type group is only a handful of lines.

    It also reports the alleles that are present but too rare to test. Those are not a failure of
    the analysis; they are the reason a patient-derived organoid biobank exists, since allele
    frequency in a cell-line panel has nothing to do with allele frequency in patients.

    Parameters
    ----------
    cancer_type : str, optional
        Cancer context, e.g. "Pancreatic Cancer" (default: "Pancreatic Cancer").
    target_mutation : str, optional
        Driver gene whose alleles define the groups (default: "KRAS").
    data_lake_path : str, optional
        Directory holding the DepMap files (default: "./data/biomni_data/data_lake").
    mutation_csv_path : str, optional
        Optional mutation table overriding the default mutation source. Must carry a ProteinChange
        column for allele resolution to work.
    min_allele_lines : int, optional
        Minimum cell lines carrying an allele before it is tested (default: 4).
    p_threshold : float, optional
        Uncorrected Welch p-value cutoff (default: 0.05).
    min_effect_difference : float, optional
        Required gene-effect difference against the comparison group; must be negative
        (default: -0.2).
    max_mutant_mean_effect : float, optional
        The allele group mean gene effect must be below this value (default: -0.3).
    max_pct_all_dependent : float, optional
        Drop genes essential in more than this percentage of all screened lines (default: 80.0).
    top_n : int, optional
        Number of allele-restricted candidates to report (default: 15).
    output_csv_path : str, optional
        Write the full allele-resolved result table to this path as well.

    Returns
    -------
    str
        A research log with the allele spectrum of the panel, a per-allele read-out of the RAS
        pathway genes the organoid study singled out, a ranked table of dependencies restricted to
        one allele and missed by the pooled test, and an explicit list of alleles that no cell-line
        panel can answer for.

    """
    import numpy as np
    import pandas as pd
    from scipy import stats

    log = [
        "=" * 78,
        f"ALLELE-RESOLVED SL DISCOVERY - {target_mutation.upper()} alleles in {cancer_type}",
        f"Run at {datetime.now(tz=_UTC).strftime('%Y-%m-%d %H:%M UTC')}",
        "=" * 78,
        "",
        f"Premise: '{target_mutation.upper()}-mutant' is a label, not a genotype - it pools biochemically",
        "different proteins, and a pooled test averages away anything restricted to one of them. The",
        "worked example is KRAS: the Nature 2026 organoid biobank found G12 organoids dependent on KRAS,",
        "EGFR and PTPN11 while Q61H organoids ignored EGFR inhibition entirely. The same argument applies",
        f"to any driver with recurrent alleles, which is why the contrast below is run on {target_mutation.upper()}.",
        "",
    ]

    try:
        bundle = _load_depmap(data_lake_path)
    except FileNotFoundError as e:
        return f"FAILURE: {e}"

    gene_effect = bundle["gene_effect"]
    symbols = [c.split(" (")[0].upper() for c in gene_effect.columns]
    cohort, match_note = _select_cancer_models(bundle["model"], cancer_type)
    cohort = cohort[cohort["ModelID"].isin(gene_effect.index)]
    if len(cohort) < 10:
        return f"FAILURE: only {len(cohort)} screened cell lines matched '{cancer_type}'; allele splitting needs more."

    try:
        cohort, mutation_source = _annotate_mutation_status(
            cohort, target_mutation, bundle["data_lake_path"], mutation_csv_path
        )
    except (RuntimeError, ValueError) as e:
        return f"FAILURE: {e}"

    cohort = cohort.copy()
    parsed = cohort.apply(
        lambda r: _normalise_allele(r["Variant"]) if r["MutationStatus"] == "MUT" else ("WT", "WT"),
        axis=1,
        result_type="expand",
    )
    cohort["allele"] = parsed[0]
    cohort["codon"] = parsed[1]
    cohort.loc[cohort["MutationStatus"] == "UNKNOWN", ["allele", "codon"]] = "UNKNOWN"

    mutant = cohort[cohort["MutationStatus"] == "MUT"]
    wildtype_ids = cohort.loc[cohort["MutationStatus"] == "WT", "ModelID"].tolist()
    allele_counts = mutant["allele"].value_counts()
    codon_counts = mutant["codon"].value_counts()

    log.append("STEP 1 | Allele spectrum of the screened panel")
    log.append(f"  Context: {match_note}")
    log.append(f"  Mutation source: {mutation_source}")
    log.append(
        f"  {len(mutant)} mutant, {len(wildtype_ids)} wild-type, "
        f"{int((cohort['MutationStatus'] == 'UNKNOWN').sum())} not profiled"
    )
    log.append("")
    log.append(f"  {'allele':<14}{'n':>4}  lines")
    log.append("  " + "-" * 70)
    for allele, count in allele_counts.items():
        names = mutant.loc[mutant["allele"] == allele, "StrippedCellLineName"].tolist()
        shown = ", ".join(names[:6]) + (f", +{len(names) - 6} more" if len(names) > 6 else "")
        log.append(f"  {allele:<14}{count:>4}  {shown}")
    if wildtype_ids:
        wt_names = cohort.loc[cohort["MutationStatus"] == "WT", "StrippedCellLineName"].tolist()
        log.append(f"  {'wild type':<14}{len(wildtype_ids):>4}  {', '.join(wt_names)}")
    log.append(f"  Codon groups: {dict(codon_counts)}")

    testable = [a for a, n in allele_counts.items() if n >= min_allele_lines and a not in {"OTHER", "MULTI"}]
    testable = [a for a in testable if not a.startswith("MULTI(")]
    untestable = [(a, int(n)) for a, n in allele_counts.items() if a not in testable]

    # --- STEP 2: the read-out the organoid study singled out -------------------------------------
    log.append("")
    log.append("STEP 2 | RAS pathway dependency, one allele at a time")
    log.append("  This is the panel the organoid study contrasted between G12 and Q61H. Values are mean")
    log.append("  gene effect; more negative means more dependent. Groups below the testing threshold are")
    log.append("  shown too, marked with *, because a descriptive mean is still worth seeing - but a mean")
    log.append("  over one or two lines is an anecdote, not an estimate.")
    log.append("")
    probe_columns = {g: gene_effect.columns[symbols.index(g)] for g in RAS_PATHWAY_PROBES if g in symbols}
    groups = [(a, mutant.loc[mutant["allele"] == a, "ModelID"].tolist()) for a, _ in allele_counts.items()]
    if wildtype_ids:
        groups.append(("WT", wildtype_ids))
    groups = [(label, ids) for label, ids in groups if len(ids) >= 1]
    header = f"  {'gene':<10}" + "".join(
        f"{(label + ('' if len(ids) >= min_allele_lines else '*')) + f'(n={len(ids)})':>14}" for label, ids in groups
    )
    log.append(header)
    log.append("  " + "-" * (len(header) - 2))
    for gene, column in probe_columns.items():
        row = f"  {gene:<10}"
        for _, ids in groups:
            values = gene_effect.loc[ids, column].dropna()
            row += f"{values.mean():>14.3f}" if len(values) else f"{'-':>14}"
        log.append(row)

    if not testable:
        log.append("")
        log.append(
            f"FAILURE: no {target_mutation.upper()} allele reaches {min_allele_lines} screened lines in this "
            "cohort, so no allele-level test can be run. The allele spectrum above is the finding: the "
            "panel cannot answer allele-level questions here."
        )
        return "\n".join(log)

    # --- STEP 3: per-allele genome-wide contrasts -----------------------------------------------
    pan_pct = (gene_effect < DEPLETION_THRESHOLD).sum() / gene_effect.notna().sum() * 100
    pan_pct.index = symbols

    def welch(group_ids, reference_ids):
        """Welch t-test of every gene between two sets of lines."""
        if len(group_ids) < 3 or len(reference_ids) < 3:
            return None
        a = gene_effect.loc[group_ids]
        b = gene_effect.loc[reference_ids]
        usable = a.columns[(a.notna().sum() >= 3) & (b.notna().sum() >= 3)]
        a, b = a[usable], b[usable]
        _, pvalue = stats.ttest_ind(a.values, b.values, axis=0, equal_var=False, nan_policy="omit")
        pvalue = np.asarray(pvalue, dtype=float)
        pvalue = np.where(np.isfinite(pvalue), pvalue, 1.0)
        group_mean = np.asarray(a.mean())
        reference_mean = np.asarray(b.mean())
        return pd.DataFrame(
            {
                "gene": [c.split(" (")[0].upper() for c in usable],
                "allele_mean": group_mean,
                "reference_mean": reference_mean,
                "difference": group_mean - reference_mean,
                "p_value": pvalue,
                "q_value": _benjamini_hochberg(pvalue),
            }
        ).set_index("gene")

    pooled = welch(mutant["ModelID"].tolist(), wildtype_ids) if len(wildtype_ids) >= 3 else None
    log.append("")
    log.append("STEP 3 | Per-allele contrasts")
    if pooled is None:
        log.append(
            f"  Pooled mutant-vs-wild-type reference: NOT AVAILABLE - only {len(wildtype_ids)} wild-type "
            "lines. Everything below is reported against the other alleles instead, which is the more "
            "informative contrast anyway."
        )
    else:
        log.append(f"  Pooled mutant-vs-wild-type reference: {len(mutant)} vs {len(wildtype_ids)} lines.")

    findings = []
    for allele in testable:
        allele_ids = mutant.loc[mutant["allele"] == allele, "ModelID"].tolist()
        other_ids = mutant.loc[mutant["allele"] != allele, "ModelID"].tolist()
        versus_other = welch(allele_ids, other_ids)
        if versus_other is None:
            log.append(f"  {allele:<8}: skipped (n={len(allele_ids)} vs {len(other_ids)} other-allele lines)")
            continue
        versus_wt = welch(allele_ids, wildtype_ids) if len(wildtype_ids) >= 3 else None
        log.append(
            f"  {allele:<8}: {len(allele_ids)} lines vs {len(other_ids)} other-allele lines"
            + (f", vs {len(wildtype_ids)} wild-type lines" if versus_wt is not None else "")
        )

        selected = versus_other[
            (versus_other["p_value"] < p_threshold)
            & (versus_other["difference"] <= min_effect_difference)
            & (versus_other["allele_mean"] <= max_mutant_mean_effect)
        ].copy()
        selected["pct_all_dependent"] = pan_pct.reindex(selected.index)
        selected = selected[selected["pct_all_dependent"] < max_pct_all_dependent]
        selected["allele"] = allele
        selected["n_allele_lines"] = len(allele_ids)
        if versus_wt is not None:
            selected = selected.join(versus_wt[["difference", "p_value"]], rsuffix="_vs_wt", how="left")
        if pooled is not None:
            selected = selected.join(pooled[["difference", "p_value"]], rsuffix="_pooled", how="left")
            # The point of the tool: keep what the pooled label could not see.
            selected["hidden_by_pooling"] = (selected["p_value_pooled"] >= p_threshold) | (
                selected["difference_pooled"] > min_effect_difference
            )
        else:
            selected["hidden_by_pooling"] = True
        findings.append(selected)

    if not findings or all(len(f) == 0 for f in findings):
        log.append("")
        log.append("No allele-restricted candidate survived the filters.")
        return "\n".join(log)

    combined = pd.concat(findings).sort_values("difference")
    hidden = combined[combined["hidden_by_pooling"]]

    log.append("")
    log.append(f"STEP 4 | Dependencies restricted to one allele (top {min(top_n, len(hidden))})")
    log.append(
        f"  {len(combined)} allele-restricted hits, of which {len(hidden)} are invisible to the pooled "
        f"{target_mutation.upper()}-mutant test."
    )
    log.append("")
    show_pooled = "p_value_pooled" in combined.columns
    header = f"  {'gene':<12}{'allele':<9}{'n':>3}{'allele eff':>12}{'other eff':>11}{'diff':>8}{'p':>10}{'q':>9}"
    if show_pooled:
        header += f"{'p(pooled)':>11}"
    header += f"{'%all dep':>10}"
    log.append(header)
    log.append("  " + "-" * (len(header) - 2))
    for gene, row in hidden.head(top_n).iterrows():
        line = (
            f"  {gene:<12}{row['allele']:<9}{int(row['n_allele_lines']):>3}{row['allele_mean']:>12.3f}"
            f"{row['reference_mean']:>11.3f}{row['difference']:>8.3f}{row['p_value']:>10.2e}{row['q_value']:>9.3f}"
        )
        if show_pooled:
            pooled_p = row.get("p_value_pooled")
            line += f"{pooled_p:>11.2e}" if pd.notna(pooled_p) else f"{'-':>11}"
        line += f"{row['pct_all_dependent']:>9.0f}%"
        log.append(line)

    genes = list(dict.fromkeys(hidden.head(top_n).index.tolist()))
    log.append("")
    log.append(f"ALLELE_RESOLVED_CANDIDATES: {', '.join(genes) if genes else '(none)'}")

    # --- STEP 5: what the panel cannot answer ---------------------------------------------------
    log.append("")
    log.append("STEP 5 | Alleles this panel cannot answer for")
    if untestable:
        for allele, count in untestable:
            names = mutant.loc[mutant["allele"] == allele, "StrippedCellLineName"].tolist()
            log.append(f"  {allele:<14}n={count:<3} ({', '.join(names)}) - below the {min_allele_lines}-line threshold")
        log.append("")
        log.append("  These are not missing data to be worked around. Cell-line panels were assembled")
        log.append("  decades ago from whatever grew on plastic, so allele frequency in the panel has no")
        log.append("  relationship to allele frequency in patients. An organoid biobank recruits from")
        log.append("  patients, so it can be built to cover a rare allele on purpose. That is the specific")
        log.append("  thing an organoid gives you here - not realism in the abstract, but the genotypes.")
    else:
        log.append(f"  None - every observed allele reached {min_allele_lines} lines.")

    if output_csv_path:
        combined.to_csv(output_csv_path)
        log.append("")
        log.append(f"Full allele-resolved table written to {output_csv_path} ({len(combined)} rows).")

    log.append("")
    log.append("HOW TO READ THIS (important)")
    log.append("  The primary contrast is allele-versus-other-alleles, not allele-versus-wild-type. It")
    log.append("  answers 'does this allele differ from the other mutants', which is the question the")
    log.append("  pooled label obscures, and it does not depend on the tiny wild-type group. It does NOT")
    log.append("  establish synthetic lethality with the driver: a gene can differ between alleles while")
    log.append("  being equally required in wild-type cells. Read the vs-wild-type columns before")
    log.append("  claiming a genotype-selective vulnerability.")
    log.append("  Subgroup testing on this many genes with this few lines will produce false positives;")
    log.append("  the q-value column is there to keep that visible. One outlying line can carry a whole")
    log.append("  allele group, so inspect the per-line values before believing any single row.")
    log.append("")
    log.append("PROVENANCE AND LIMITS")
    for item in bundle["provenance"]:
        log.append(f"  - {item}")
    log.append(f"  - Mutation calls: {mutation_source}")
    log.append(
        "  - Allele assignment uses the protein change string only. Zygosity, allelic imbalance and "
        "expression of the mutant allele are not modelled, and a line carrying two distinct changes is "
        "excluded rather than assigned."
    )
    log.append(
        "  - No organoid data is involved: 0 of the 24 organoid models in the local DepMap snapshot have "
        "CRISPR or expression data."
    )
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
    serum_contrast = _serum_free_contrast(data_lake_path)
    serum_programme = {}
    for programme_name, programme_genes in SERUM_DEPENDENT_PROGRAMS.items():
        for programme_gene in programme_genes:
            serum_programme.setdefault(programme_gene, programme_name)

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
    if serum_contrast is None:
        log.append("No serum-free lines annotated - the serum-masking check is skipped.")
    else:
        log.append(
            f"Serum axis reference: {serum_contrast['n_serum_free']} serum-free lines, each z-scored "
            "against serum-cultured lines of its own lineage."
        )

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
                masked_by.append((factor, spec["note"], spec.get("withholds", False)))
            if gene in spec["preserved_downstream"]:
                preserved_by.append(factor)
        if masked_by:
            for factor, note, _ in masked_by:
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

        # --- serum masking (measured) ---
        # Organoid medium is serum-free, so anything the cell normally scavenges from FBS has to be
        # made instead. This is measured the same way anchorage is, on the lines DepMap annotates.
        serum_z = float("nan")
        if serum_contrast is not None and gene in serum_contrast["z"].index:
            serum_z = float(serum_contrast["z"].loc[gene])
            scored = int(serum_contrast["n_models_scored"].get(gene, 0))
            programme = serum_programme.get(gene)
            if serum_z <= -1.0:
                reading = "stronger without serum - organoid medium should unmask it"
            elif serum_z >= 1.0:
                reading = "weaker without serum - organoid medium should blunt it"
            else:
                reading = "no meaningful shift"
            log.append(
                f"  Serum-free sensitivity: lineage-matched z={serum_z:+.2f} over {scored} serum-free "
                f"lines - {reading}" + (f" [{programme}]" if programme else "")
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
            factor, _, withholds = masked_by[0]
            if withholds:
                verdict = (
                    f"ORGANOID-MASKED - the dependency needs something organoid medium does not contain "
                    f"({factor}), so there is nothing for this gene to act on and the knockout will look "
                    "inert. To ask the question at all, add the missing component back (lipoprotein-"
                    "supplemented or serum-containing medium); otherwise the organoid returns a false negative."
                )
            else:
                verdict = (
                    f"ORGANOID-MASKED - standard medium supplies {factor}, which rescues this dependency. "
                    f"Drop {factor} from the medium (or use a reduced-factor formulation) or the organoid "
                    "will return a false negative."
                )
        elif serum_z <= -1.0:
            verdict = (
                f"ORGANOID-ENHANCED (lipid axis) - measurably stronger in serum-free culture "
                f"(z={serum_z:+.2f}), which is how organoids are fed. The 2D screen understated this "
                "dependency because the serum supplied what the cell would otherwise have to make. "
                "Cheapest confirmation is lipid-depleted serum in the same cell lines, before any organoid."
            )
        elif anchorage_p < 1e-3 and anchorage_difference < -0.15:
            verdict = (
                "ORGANOID-ENHANCED - the dependency is measurably stronger in anchored culture, so a "
                "3D matrix-embedded organoid should show a larger effect than the 2D screen did."
            )
        elif serum_z >= 1.0:
            verdict = (
                f"ORGANOID-WEAKENED (lipid axis) - weaker in serum-free culture (z={serum_z:+.2f}). "
                "Expect the organoid to understate this dependency relative to the 2D screen."
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
    if serum_contrast is not None:
        log.append(
            f"  - Serum-free sensitivity is measured on only {serum_contrast['n_serum_free']} lines, none "
            "of them organoids and none from this lineage; per-gene z-scores are noisy. "
            "discover_serum_masked_dependencies() reports the better-powered programme-level test."
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


# --- rank_organoid_sl_candidates decision rules -------------------------------------------------
# Every threshold below is a stated rule, not a fitted parameter, and each is module-level so a
# reader can see - and change - the exact line at which a candidate is reclassified.
#
# The rule that matters most is SL_GENOTYPE_DIFFERENCE. The three discovery axes each answer their
# own question (does a subtype hide it, does the allele label dilute it, does serum mask it) and
# none of them asks whether the dependency is genotype-selective. discover_serum_masked_dependencies
# in particular computes the mutant-minus-wild-type difference and prints it as a column, but sorts
# and filters purely on the serum z-score - so a gene required by every cell regardless of driver
# status reaches the candidate list on equal footing with a genuine synthetic lethal. That is the
# organoid-specific core fitness class the Nature 2026 biobank reported, and it is not a drug
# target: normal tissue in a patient is not lipid-starved, so the matched normal organoid will die
# alongside the tumour organoid. This gate separates the two.
SL_GENOTYPE_DIFFERENCE = -0.15  # mutant minus wild-type gene effect required to call it selective
SL_WINDOW_MAX_PCT_ALL = 50.0  # above this share of all screened lines the gene is drifting pan-essential
SL_WINDOW_WATCH_PCT_ALL = 35.0  # softer band, penalised but not disqualifying
SL_WINDOW_MAX_NORMAL_TPM = 20.0  # GTEx normal-tissue expression above which the window must be proven
PARALOG_CONFOUND_R = 0.25  # identical to the PARALOG ALERT threshold in check_dependency_confounders
PARALOG_WATCH_R = 0.15  # below the alert, but high enough that the rescue arm must control for it
ANCHORAGE_INSTRUMENT_DIFF = -0.15  # matches the ORGANOID-ENHANCED rule in assess_organoid_transferability
ANCHORAGE_INSTRUMENT_P = 1e-3
SERUM_UNMASK_Z = -1.0  # matches the lipid-axis rule in assess_organoid_transferability


def _paralog_correlation(gene: str, gene_effect, gene_columns, expression_bundle):
    """Strongest correlation between a paralog's expression and this gene's dependency.

    Reuses the family-prefix definition from ``check_dependency_confounders`` so the two tools
    cannot disagree about what counts as a paralog. Only the family is correlated rather than the
    whole transcriptome, because the ranking needs the family maximum, not a genome-wide rank.
    """
    if expression_bundle is None or gene not in gene_columns:
        return None, float("nan")
    expression = expression_bundle["expression"]
    expression_columns = expression_bundle["gene_columns"]

    family_prefixes = {re.sub(r"(?<=\d)[A-Z]$", "", gene), re.sub(r"\d+$", "", gene)}
    family_prefixes = {p for p in family_prefixes if len(p) >= 3 and p != gene}
    family = sorted(
        {
            symbol
            for symbol in expression_columns
            if symbol != gene and any(symbol.startswith(prefix) for prefix in family_prefixes)
        }
    )
    if not family:
        return None, float("nan")

    dependency = gene_effect[gene_columns[gene]].dropna()
    shared = dependency.index.intersection(expression.index)
    if len(shared) < 30:
        return None, float("nan")
    dependency = dependency.loc[shared]

    best_symbol, best_r = None, 0.0
    for symbol in family:
        values = expression.loc[shared, expression_columns[symbol]]
        if values.notna().sum() < 30 or values.std() == 0:
            continue
        r = float(dependency.corr(values))
        if r == r and abs(r) > abs(best_r):
            best_symbol, best_r = symbol, r
    return best_symbol, best_r if best_symbol else float("nan")


def rank_organoid_sl_candidates(
    candidate_genes,
    target_mutation: str = "KRAS",
    cancer_type: str = "Pancreatic Cancer",
    axis_origins: dict | None = None,
    literature_scores: dict | None = None,
    q_values: dict | None = None,
    normal_tissue: str = "Pancreas",
    data_lake_path: str | None = None,
    mutation_csv_path: str | None = None,
    top_n: int = 25,
) -> str:
    """Rank candidates by how well an organoid experiment could establish synthetic lethality.

    The other tools in this module each report one axis and never compare candidates with one
    another: ``assess_organoid_transferability`` emits an independent verdict per gene, and the
    three discovery tools emit independent candidate lists. Deciding which candidate an organoid
    should actually be spent on therefore had to be done by hand, outside the tooling. This function
    does it inside, from explicit rules, so the conclusion is reproducible and auditable rather than
    a matter of who read the reports.

    Two distinctions drive the ranking, and neither is made by the per-axis tools:

    1. **Genotype-selective versus core fitness.** A dependency that is equally strong in
       driver-wild-type cells is not synthetic lethal, however dramatically an organoid unmasks it.
       The serum axis is where this matters: it selects on the serum-free z-score alone, so
       organoid-specific core fitness genes arrive looking like discoveries. Candidates that fail
       ``SL_GENOTYPE_DIFFERENCE`` are reclassified CORE-FITNESS-NOT-SL and cannot be top-ranked.
    2. **Why the organoid helps.** ``ORGANOID-ENHANCED`` is awarded for two unrelated reasons -
       anchorage in matrix, or the absence of serum. Only the first means the organoid format is a
       better instrument for the same question. The second means the medium changes which metabolic
       genes are limiting, which in a fed patient is a window risk rather than a target.

    Parameters
    ----------
    candidate_genes : list[str] | str
        Candidates to rank; a list, or a comma/whitespace separated string.
    target_mutation : str, optional
        The driver gene defining the genotype contrast (default: "KRAS").
    cancer_type : str, optional
        Cancer context used to select the cell-line cohort (default: "Pancreatic Cancer").
    axis_origins : dict, optional
        Mapping of gene to the discovery axes it came from, e.g. ``{"SCAP": ["allele", "lipid"]}``.
        Convergence across independent axes earns a bounded bonus.
    literature_scores : dict, optional
        Mapping of gene to the 0-100 score from ``validate_sl_candidates_with_pubmed``. Used only as
        a modulator; absence is treated as neutral, never as evidence against.
    q_values : dict, optional
        Mapping of gene to the genome-wide FDR q-value from ``discover_synthetic_lethal_candidates``.
        Omitted rather than recomputed, because an FDR over a handful of candidates would not mean
        what a genome-wide one does.
    normal_tissue : str, optional
        GTEx tissue used as the normal counterpart for the window check (default: "Pancreas").
    data_lake_path : str, optional
        Directory holding the DepMap and GTEx files.
    mutation_csv_path : str, optional
        Optional mutation table overriding the default mutation source.
    top_n : int, optional
        Number of ranked rows to print (default: 25).

    Returns
    -------
    str
        A research log with the rules applied, a per-candidate evidence block, a ranked table with a
        class and an organoid reading for every candidate, and machine-readable
        ``ORGANOID_SL_RANKING`` / ``ORGANOID_TOP_PICK`` lines.

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

    cohort, match_note = _select_cancer_models(bundle["model"], cancer_type)
    cohort = cohort[cohort["ModelID"].isin(gene_effect.index)]
    if len(cohort) < 5:
        return f"FAILURE: only {len(cohort)} screened cell lines matched '{cancer_type}'."
    try:
        cohort, mutation_source = _annotate_mutation_status(
            cohort, target_mutation, bundle["data_lake_path"], mutation_csv_path
        )
    except (RuntimeError, ValueError) as e:
        return f"FAILURE: {e}"
    mutant_ids = cohort.loc[cohort["MutationStatus"] == "MUT", "ModelID"]
    wildtype_ids = cohort.loc[cohort["MutationStatus"] == "WT", "ModelID"]

    try:
        expression_bundle = _load_expression(data_lake_path)
    except FileNotFoundError:
        expression_bundle = None
    gtex = _load_gtex(data_lake_path)
    serum_contrast = _serum_free_contrast(data_lake_path)
    axis_origins = {k.upper(): v for k, v in (axis_origins or {}).items()}
    literature_scores = {k.upper(): v for k, v in (literature_scores or {}).items()}
    q_values = {k.upper(): v for k, v in (q_values or {}).items()}

    driver = target_mutation.strip().upper()
    log = [
        "=" * 78,
        f"ORGANOID SL RANKING - {driver}-mutant {cancer_type}",
        f"Run at {datetime.now(tz=_UTC).strftime('%Y-%m-%d %H:%M UTC')}",
        "=" * 78,
        "",
        "Question: of these candidates, which one could a patient-derived organoid actually",
        "establish as synthetic lethal - not merely reproduce, and not confuse with core fitness.",
        "",
        "RULES APPLIED (all deterministic; change the constants at the top of this section to change",
        "the ranking - no judgement is exercised anywhere below)",
        f"  Gate 1  genotype selectivity : mutant-minus-wild-type gene effect <= {SL_GENOTYPE_DIFFERENCE}",
        "          A gene that fails this is required regardless of driver status. An organoid will",
        "          reproduce it faithfully AND kill the matched normal organoid. Class:",
        "          CORE-FITNESS-NOT-SL. It cannot be top-ranked whatever else it scores.",
        f"  Gate 2  therapeutic window   : < {SL_WINDOW_MAX_PCT_ALL:.0f}% of all screened lines dependent,",
        f"          < {SL_WINDOW_MAX_NORMAL_TPM:.0f} TPM in normal {normal_tissue}, 0% of wild-type lines dependent",
        f"  Gate 3  paralog robustness   : family-max |r| >= {PARALOG_CONFOUND_R} confounds the call;",
        f"          >= {PARALOG_WATCH_R} is a watch flag that the rescue arm must control for",
        "  Signal  organoid reading     : ORGANOID-INSTRUMENT (anchorage or niche coupling - the 3D",
        "          format measures the same question better) vs ORGANOID-MEDIUM-EFFECT (serum-free",
        "          unmasking - the medium changes the question) vs ORGANOID-NEUTRAL (restates 2D)",
        "  Bonus   axis convergence, bounded; literature score, bounded and never negative",
        "",
        f"Cohort: {match_note}; {len(mutant_ids)} {driver}-mutant, {len(wildtype_ids)} wild-type screened lines.",
        f"Mutation source: {mutation_source}",
        f"Anchorage reference panel: {len(adherent)} adherent vs {len(suspension)} suspension lines.",
    ]
    if serum_contrast is None:
        log.append("Serum axis unavailable - the medium reading is skipped.")
    if gtex is None:
        log.append(f"GTEx unavailable - the normal {normal_tissue} window check is skipped.")
    if expression_bundle is None:
        log.append("Expression matrix unavailable - the paralog check is skipped.")

    rows = []
    for gene in genes:
        if gene not in columns:
            log.append("")
            log.append(f"### {gene}")
            log.append("  Not in the CRISPR library; skipped.")
            continue
        series = gene_effect[columns[gene]]

        mutant = series.reindex(mutant_ids).dropna()
        wildtype = series.reindex(wildtype_ids).dropna()
        mutant_mean = float(mutant.mean()) if len(mutant) else float("nan")
        wildtype_mean = float(wildtype.mean()) if len(wildtype) else float("nan")
        difference = mutant_mean - wildtype_mean if len(mutant) and len(wildtype) else float("nan")
        genotype_p = float("nan")
        if len(mutant) >= 3 and len(wildtype) >= 3:
            genotype_p = float(stats.ttest_ind(mutant, wildtype, equal_var=False)[1])
        pct_mutant = 100.0 * (mutant <= DEPLETION_THRESHOLD).mean() if len(mutant) else float("nan")
        pct_wildtype = 100.0 * (wildtype <= DEPLETION_THRESHOLD).mean() if len(wildtype) else float("nan")
        pct_all = 100.0 * (series.dropna() <= DEPLETION_THRESHOLD).mean()

        a = series.reindex(adherent).dropna()
        s = series.reindex(suspension).dropna()
        anchorage_difference, anchorage_p = float("nan"), 1.0
        if len(a) > 10 and len(s) > 10:
            anchorage_difference = float(a.mean() - s.mean())
            anchorage_p = float(stats.ttest_ind(a, s, equal_var=False)[1])

        serum_z = float("nan")
        if serum_contrast is not None and gene in serum_contrast["z"].index:
            serum_z = float(serum_contrast["z"].loc[gene])

        niche_masked = [f for f, spec in NICHE_FACTORS.items() if gene in spec["masked_upstream"]]
        niche_preserved = [f for f, spec in NICHE_FACTORS.items() if gene in spec["preserved_downstream"]]

        normal_tpm = float("nan")
        if gtex is not None:
            hit = gtex[(gtex.Tissue == normal_tissue) & (gtex.Gene == gene)]
            if len(hit):
                normal_tpm = float(hit.Expression.iloc[0])

        paralog_symbol, paralog_r = _paralog_correlation(gene, gene_effect, columns, expression_bundle)

        # --- classification ---------------------------------------------------------------------
        contradictions = []
        genotype_selective = difference == difference and difference <= SL_GENOTYPE_DIFFERENCE
        if not genotype_selective:
            shown = f"{difference:+.3f}" if difference == difference else "not measurable"
            contradictions.append(
                f"mutant-minus-wild-type {shown} does not reach {SL_GENOTYPE_DIFFERENCE} - the "
                f"dependency does not track {driver} status"
            )

        paralog_confounded = paralog_r == paralog_r and abs(paralog_r) >= PARALOG_CONFOUND_R
        paralog_watch = paralog_r == paralog_r and PARALOG_WATCH_R <= abs(paralog_r) < PARALOG_CONFOUND_R
        if paralog_confounded:
            contradictions.append(
                f"paralog {paralog_symbol} expression tracks the dependency (r={paralog_r:+.3f}) - "
                "this reads as a paralog-loss effect, not a driver interaction"
            )
        elif paralog_watch:
            contradictions.append(
                f"paralog {paralog_symbol} r={paralog_r:+.3f} sits below the alert threshold but above "
                "background; stratify by its expression and keep it in the rescue arm"
            )

        window_failed = False
        if pct_all >= SL_WINDOW_MAX_PCT_ALL:
            window_failed = True
            contradictions.append(f"{pct_all:.0f}% of all screened lines dependent - no therapeutic window")
        elif pct_all >= SL_WINDOW_WATCH_PCT_ALL:
            contradictions.append(f"{pct_all:.0f}% of all screened lines dependent - window is narrow")
        if normal_tpm == normal_tpm and normal_tpm >= SL_WINDOW_MAX_NORMAL_TPM:
            contradictions.append(
                f"{normal_tpm:.0f} TPM in normal {normal_tissue} - the matched normal organoid arm must "
                "run before any efficacy claim"
            )
        if pct_wildtype == pct_wildtype and pct_wildtype > 0:
            contradictions.append(f"{pct_wildtype:.0f}% of wild-type lines also dependent - not fully selective")
        q_value = q_values.get(gene)
        if q_value is not None and q_value > 0.1:
            contradictions.append(f"genome-wide FDR q={q_value:.3f} - weak multiple-testing support")

        if not genotype_selective:
            candidate_class = "CORE-FITNESS-NOT-SL"
        elif paralog_confounded:
            candidate_class = "PARALOG-CONFOUNDED"
        elif window_failed:
            candidate_class = "NO-WINDOW"
        else:
            candidate_class = "SL-CANDIDATE"

        # --- organoid reading -------------------------------------------------------------------
        anchorage_instrument = (
            anchorage_p < ANCHORAGE_INSTRUMENT_P
            and anchorage_difference == anchorage_difference
            and anchorage_difference <= ANCHORAGE_INSTRUMENT_DIFF
        )
        if niche_masked:
            organoid_reading = "ORGANOID-MASKED"
            organoid_note = (
                f"standard medium supplies or withholds {niche_masked[0]}, which changes the readout; "
                "modify the medium or the organoid returns a false negative"
            )
        elif anchorage_instrument or niche_preserved:
            organoid_reading = "ORGANOID-INSTRUMENT"
            reason = (
                f"anchorage diff {anchorage_difference:+.3f}, p={anchorage_p:.1e}"
                if anchorage_instrument
                else f"downstream of {', '.join(niche_preserved)}, not rescued by the medium"
            )
            organoid_note = f"the 3D format amplifies the same dependency ({reason}) - an organoid measures it better"
        elif serum_z == serum_z and serum_z <= SERUM_UNMASK_Z:
            organoid_reading = "ORGANOID-MEDIUM-EFFECT"
            organoid_note = (
                f"serum-free z={serum_z:+.2f}: the organoid medium changes which metabolic genes are "
                "limiting. Falsify with lipid-depleted serum in the same cell lines first - it is far "
                "cheaper, and a patient's normal tissue is not lipid-starved"
            )
        else:
            organoid_reading = "ORGANOID-NEUTRAL"
            organoid_note = (
                "no medium conflict and no anchorage effect - the organoid should restate the "
                "cell-line result, so it buys subtype and patient diversity, not mechanism"
            )

        # --- score ------------------------------------------------------------------------------
        # Directness-weighted, in the same spirit as generate_sl_evidence_dossier: the genotype
        # contrast dominates, the organoid reading modulates, everything else can only trim.
        score = 0.0
        if genotype_selective:
            score += min(40.0, abs(difference) * 100.0)
            score += min(15.0, max(0.0, (pct_mutant - pct_all)) * 0.5) if pct_mutant == pct_mutant else 0.0
        score += {"ORGANOID-INSTRUMENT": 25.0, "ORGANOID-NEUTRAL": 8.0, "ORGANOID-MASKED": 5.0}.get(
            organoid_reading, 0.0
        )
        axes = axis_origins.get(gene, [])
        score += min(10.0, max(0, len(axes) - 1) * 5.0)
        literature = literature_scores.get(gene)
        if literature is not None:
            score += min(15.0, literature * 0.15)
        if paralog_confounded:
            score -= 30.0
        elif paralog_watch:
            score -= 12.0
        if pct_all >= SL_WINDOW_MAX_PCT_ALL:
            score -= 20.0
        elif pct_all >= SL_WINDOW_WATCH_PCT_ALL:
            score -= 8.0
        if normal_tpm == normal_tpm and normal_tpm >= SL_WINDOW_MAX_NORMAL_TPM:
            score -= 10.0
        if pct_wildtype == pct_wildtype and pct_wildtype > 0:
            score -= 10.0
        if q_value is not None and q_value > 0.1:
            score -= 5.0
        score = max(0.0, score)

        rows.append(
            {
                "gene": gene,
                "class": candidate_class,
                "organoid_reading": organoid_reading,
                "organoid_note": organoid_note,
                "score": score,
                "difference": difference,
                "genotype_p": genotype_p,
                "pct_mutant": pct_mutant,
                "pct_wildtype": pct_wildtype,
                "pct_all": pct_all,
                "anchorage_difference": anchorage_difference,
                "anchorage_p": anchorage_p,
                "serum_z": serum_z,
                "normal_tpm": normal_tpm,
                "paralog_symbol": paralog_symbol,
                "paralog_r": paralog_r,
                "axes": axes,
                "literature": literature,
                "q_value": q_value,
                "contradictions": contradictions,
            }
        )

    if not rows:
        return "FAILURE: none of the supplied genes are present in the CRISPR library."

    # SL-CANDIDATE outranks everything else regardless of score, so that a spectacular core fitness
    # gene can never displace a modest but genuine synthetic lethal.
    class_order = {"SL-CANDIDATE": 0, "NO-WINDOW": 1, "PARALOG-CONFOUNDED": 2, "CORE-FITNESS-NOT-SL": 3}
    rows.sort(key=lambda r: (class_order[r["class"]], -r["score"]))

    log.append("")
    log.append("=" * 78)
    log.append("RANKING")
    log.append("=" * 78)
    header = (
        f"  {'gene':<11}{'score':>6}{'MUT-WT':>8}{'%MUTdep':>8}{'%WTdep':>7}{'%ALLdep':>8}"
        f"{'anchor':>8}{'serumZ':>7}{'nTPM':>7}  {'class':<20}organoid reading"
    )
    log.append(header)
    log.append("  " + "-" * (len(header) + 6))
    for row in rows[:top_n]:

        def show(value, fmt):
            """Format a cell, keeping the column width when the value is missing."""
            if value == value:
                return format(value, fmt)
            width = int(re.search(r"(\d+)", fmt).group(1))
            return "-".rjust(width)

        log.append(
            f"  {row['gene']:<11}{row['score']:>6.1f}{show(row['difference'], '>+8.3f')}"
            f"{show(row['pct_mutant'], '>7.0f')}%{show(row['pct_wildtype'], '>6.0f')}%"
            f"{show(row['pct_all'], '>7.0f')}%{show(row['anchorage_difference'], '>8.3f')}"
            f"{show(row['serum_z'], '>7.2f')}{show(row['normal_tpm'], '>7.1f')}  "
            f"{row['class']:<20}{row['organoid_reading']}"
        )

    log.append("")
    log.append("=" * 78)
    log.append("PER-CANDIDATE READING")
    log.append("=" * 78)
    for row in rows[:top_n]:
        log.append("")
        log.append(f"### {row['gene']}  [{row['class']}]  score {row['score']:.1f}")
        if row["axes"]:
            log.append(f"  Discovery axes: {', '.join(row['axes'])}")
        if row["literature"] is not None:
            log.append(f"  Literature score: {row['literature']}/100")
        log.append(f"  Organoid reading: {row['organoid_reading']} - {row['organoid_note']}")
        if row["contradictions"]:
            for item in row["contradictions"]:
                log.append(f"  Against: {item}")
        else:
            log.append("  Against: nothing in these checks argues against it.")
        if row["class"] == "SL-CANDIDATE" and row["organoid_reading"] == "ORGANOID-INSTRUMENT":
            verdict = (
                "VERDICT: TAKE TO ORGANOID - genotype-selective with a window, and the 3D format "
                "measures it better than the 2D screen did. This is what an organoid is for."
            )
        elif row["class"] == "SL-CANDIDATE" and row["organoid_reading"] == "ORGANOID-MEDIUM-EFFECT":
            verdict = (
                "VERDICT: FALSIFY IN 2D FIRST - genotype-selective, but the organoid signal would come "
                "from the medium. Lipid-depleted serum in the same cell lines answers it more cheaply."
            )
        elif row["class"] == "SL-CANDIDATE":
            verdict = (
                "VERDICT: ORGANOID ADDS LITTLE - genotype-selective, but the organoid would restate the "
                "cell-line result. Use it for subtype and patient diversity, not for mechanism."
            )
        elif row["class"] == "CORE-FITNESS-NOT-SL":
            verdict = (
                f"VERDICT: NOT A SYNTHETIC LETHAL - required irrespective of {driver} status. An organoid "
                "would show a large effect and so would the matched normal organoid."
            )
        elif row["class"] == "PARALOG-CONFOUNDED":
            verdict = (
                "VERDICT: STRATIFY FIRST - the paralog explains the dependency better than the driver. "
                "Re-run the contrast within paralog-high and paralog-low lines before any organoid."
            )
        else:
            verdict = (
                "VERDICT: NO WINDOW - too many unselected lines depend on it for a tumour-selective "
                "effect to be plausible."
            )
        log.append(f"  {verdict}")

    ranked_sl = [row["gene"] for row in rows if row["class"] == "SL-CANDIDATE"]
    take_to_organoid = [
        row["gene"]
        for row in rows
        if row["class"] == "SL-CANDIDATE" and row["organoid_reading"] == "ORGANOID-INSTRUMENT"
    ]
    demoted = [row["gene"] for row in rows if row["class"] == "CORE-FITNESS-NOT-SL"]

    log.append("")
    log.append("=" * 78)
    log.append("CONCLUSION")
    log.append("=" * 78)
    if take_to_organoid:
        log.append(f"  Spend the organoid on: {take_to_organoid[0]}")
        log.append(
            "  It is the highest-scoring candidate that is both genotype-selective and measurably"
        )
        log.append("  stronger in anchored 3D culture, which is the only combination where the organoid")
        log.append("  answers a question the cell-line panel cannot.")
        if len(take_to_organoid) > 1:
            log.append(f"  Same category, as mechanism controls: {', '.join(take_to_organoid[1:])}")
    elif ranked_sl:
        log.append(f"  No candidate is both genotype-selective and organoid-enhanced. Best available: {ranked_sl[0]}")
        log.append("  An organoid would restate the cell-line result for it; spend on patient diversity, not mechanism.")
    else:
        log.append("  No candidate passes the genotype-selectivity gate. Nothing here is a synthetic lethal.")
    if demoted:
        log.append("")
        log.append(f"  Demoted as core fitness rather than synthetic lethal ({len(demoted)}):")
        log.append(f"    {', '.join(demoted)}")
        log.append("    These would score in an organoid, and in a normal organoid too. The serum axis")
        log.append("    selects on the medium effect alone, which is why they reached the candidate list.")

    log.append("")
    log.append(f"ORGANOID_SL_RANKING: {', '.join(ranked_sl) if ranked_sl else '(none)'}")
    log.append(f"ORGANOID_TOP_PICK: {take_to_organoid[0] if take_to_organoid else (ranked_sl[0] if ranked_sl else '(none)')}")
    log.append(f"ORGANOID_CORE_FITNESS_DEMOTED: {', '.join(demoted) if demoted else '(none)'}")
    # One pipe-delimited row per candidate, in ranked order, so a caller can build its own summary
    # table without re-parsing the human-readable RANKING block above.
    # gene|score|class|organoid reading|mutant-minus-wild-type|%all dependent|normal-tissue TPM
    for row in rows:
        def field(value):
            return f"{value:.3f}" if value == value else ""

        log.append(
            f"ORGANOID_SL_ROW: {row['gene']}|{row['score']:.1f}|{row['class']}|{row['organoid_reading']}"
            f"|{field(row['difference'])}|{field(row['pct_all'])}|{field(row['normal_tpm'])}"
        )

    log.append("")
    log.append("PROVENANCE AND LIMITS")
    for item in bundle["provenance"]:
        log.append(f"  - {item}")
    log.append(f"  - Mutation calls: {mutation_source}")
    if expression_bundle is not None:
        log.append(f"  - {expression_bundle['provenance']}")
    log.append(
        "  - Every number here is measured on cell lines. No organoid dependency is measured "
        "anywhere in this module: 0 of the 24 organoid models in the local DepMap snapshot carry "
        "CRISPR or expression data."
    )
    log.append(
        "  - The anchorage and serum axes are proxies for two properties of organoid culture. They "
        "are not the culture, and a proxy that ranks candidates well can still be wrong about any "
        "individual one."
    )
    log.append(
        "  - The score is an ordering device, not a probability. The class and the contradictions "
        "are the parts to read; two candidates a few points apart are not distinguishable."
    )
    return "\n".join(log)
