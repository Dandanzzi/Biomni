"""Translational tool chain for pancreatic cancer (PDAC): driver landscape -> drug -> organoid -> verdict.

The synthetic lethality modules in this package stop where a wet lab starts. They take a driver the
user already named, rank candidate partners from DepMap, and (in ``organoid_sl``) write a protocol.
Nothing reads back what the experiment measured, so the one question a translational project exists
to answer - *was the prediction right?* - was left to be settled by eye, outside the tooling.

This module closes that loop for the PDAC workflow:

1. ``profile_pdac_driver_landscape``  - which drivers are actually testable in this cohort, before
   any discovery run commits to one. Recurrence alone is useless: KRAS is mutated in ~90% of PDAC
   lines, which is exactly why a KRAS-wild-type control group barely exists.
2. ``map_sl_candidates_to_drugs``     - candidate gene -> inhibitor, with clinical phase and MOA, and
   an explicit "no chemical probe, CRISPR-only" verdict for the genes that have none.
3. ``analyze_organoid_drug_response`` - fits measured PDO dose-response curves (4-parameter logistic),
   reports IC50 / AUC / Emax per organoid and compares them between genotype groups.
4. ``analyze_crispr_validation``      - turns measured organoid CRISPR knockout viability into a
   per-gene effect with a genotype-selectivity test.
5. ``compare_prediction_with_experiment`` - scores the predictions against the measurements:
   rank concordance, precision@k, and a named list of which predictions were confirmed and which
   were falsified.

Design notes
------------
* Stages 3-5 read CSV files the lab produces. Column names are parameters, so the tools fit an
  existing LIMS export rather than forcing a schema; the defaults are documented per function.
* No stage ever reports a prediction as validated on the strength of a direction alone - an effect
  size, a p-value and the group sizes are always printed next to the verdict.
* Negative results are first-class output. ``compare_prediction_with_experiment`` prints falsified
  predictions before confirmed ones, because that list is what a prediction method is judged on.
"""

import os
from datetime import datetime

from biomni.tool.synthetic_lethality import (
    _UTC,
    _annotate_copy_number_status,
    _annotate_mutation_status,
    _load_depmap,
    _parse_gene_list,
    _resolve_data_lake,
    _select_cancer_models,
    _synlethdb_partners,
)

# Canonical PDAC drivers. Sources: TCGA PAAD and ICGC PACA consensus driver lists; the point of
# hard-coding them is that a cohort of ~60 cell lines cannot rediscover drivers de novo, so the
# tool tests a curated list for *testability* rather than pretending to nominate drivers.
PDAC_DRIVERS_MUTATION = [
    "KRAS", "TP53", "CDKN2A", "SMAD4", "ARID1A", "RNF43", "GNAS", "TGFBR2", "ACVR1B",
    "KDM6A", "KMT2C", "KMT2D", "BRCA1", "BRCA2", "PALB2", "ATM", "STK11", "RBM10",
]
PDAC_DRIVERS_AMPLIFICATION = ["MYC", "GATA6", "ERBB2", "CCNE1", "AKT2"]
PDAC_DRIVERS_DELETION = ["CDKN2A", "SMAD4", "TP53"]

# A Welch comparison needs this many lines per group before its p-value means anything here.
MIN_GROUP_SIZE = 3
INFORMATIVE_GROUP_SIZE = 8

CLINICAL_PHASE_RANK = {
    "Launched": 5,
    "Phase 3": 4,
    "Phase 2/Phase 3": 4,
    "Phase 2": 3,
    "Phase 1/Phase 2": 2,
    "Phase 1": 2,
    "Preclinical": 1,
    "Withdrawn": 0,
}

# A Repurposing Hub "target" annotation does not mean the compound inhibits the protein - cystine
# is annotated to SLC7A11 (its transporter) and glutathione to GPX4 (its substrate). Only an
# inhibitory mechanism makes a compound a usable stand-in for the knockout in an organoid arm.
_INHIBITORY_MOA_TERMS = ("inhibitor", "antagonist", "blocker", "degrader", "inverse agonist", "disruptor")

_DRUG_TABLE_CACHE: dict = {}


def _is_inhibitory(moa: str) -> bool:
    text = str(moa).lower()
    return any(term in text for term in _INHIBITORY_MOA_TERMS)


def _load_repurposing_hub(data_lake_path: str | None = None):
    """Load the Broad Repurposing Hub table (pert_iname, clinical_phase, moa, target, indication)."""
    import pandas as pd

    resolved = _resolve_data_lake(data_lake_path)
    if resolved in _DRUG_TABLE_CACHE:
        return _DRUG_TABLE_CACHE[resolved]

    path = os.path.join(resolved, "broad_repurposing_hub_phase_moa_target_info.parquet")
    if not os.path.exists(path):
        _DRUG_TABLE_CACHE[resolved] = None
        return None

    table = pd.read_parquet(path)
    records = []
    for _, row in table.iterrows():
        targets = str(row.get("target") or "")
        if not targets or targets == "nan" or targets == "None":
            continue
        phase = str(row.get("clinical_phase") or "Preclinical").strip() or "Preclinical"
        for symbol in targets.split("|"):
            symbol = symbol.strip().upper()
            if symbol:
                records.append(
                    {
                        "gene": symbol,
                        "drug": str(row.get("pert_iname") or ""),
                        "clinical_phase": phase,
                        "phase_rank": CLINICAL_PHASE_RANK.get(phase, 1),
                        "moa": str(row.get("moa") or ""),
                        "indication": str(row.get("indication") or ""),
                        "disease_area": str(row.get("disease_area") or ""),
                    }
                )
    frame = pd.DataFrame(records)
    bundle = {"table": frame, "path": path, "n_drugs": table.shape[0], "n_genes": frame["gene"].nunique()}
    _DRUG_TABLE_CACHE[resolved] = bundle
    return bundle


def _read_table(path: str, required: dict, label: str):
    """Read a CSV/TSV and verify the caller's column mapping. Raises ValueError with what is missing."""
    import pandas as pd

    if not os.path.exists(path):
        raise ValueError(f"{label} file not found: {path}")
    separator = "\t" if path.lower().endswith((".tsv", ".tab")) else ","
    frame = pd.read_csv(path, sep=separator)
    missing = {role: column for role, column in required.items() if column and column not in frame.columns}
    if missing:
        raise ValueError(
            f"{label} ({path}) is missing column(s) {sorted(missing.values())} for {sorted(missing)}. "
            f"Columns present: {list(frame.columns)}. Pass the matching *_column arguments."
        )
    return frame


# ---------------------------------------------------------------------------
# Stage 1: which PDAC driver can this cohort actually test?
# ---------------------------------------------------------------------------
def profile_pdac_driver_landscape(
    cancer_type: str = "Pancreatic Cancer",
    mutation_drivers=None,
    amplification_drivers=None,
    deletion_drivers=None,
    data_lake_path: str | None = None,
    mutation_csv_path: str | None = None,
    output_csv_path: str | None = None,
) -> str:
    """Report which PDAC driver alterations split the DepMap cohort into testable groups.

    Every downstream synthetic-lethality tool needs a driver and an implicit assumption that both a
    mutant and a wild-type group exist. In PDAC that assumption fails for the most famous driver:
    KRAS is altered in the overwhelming majority of lines, so the "wild-type" arm is a handful of
    atypical models. This tool makes that visible before a discovery run is committed to, by
    counting the groups each candidate driver would produce and labelling the comparison TESTABLE,
    UNDERPOWERED or UNTESTABLE.

    Parameters
    ----------
    cancer_type : str, optional
        Cancer context passed to the DepMap lineage matcher (default: "Pancreatic Cancer").
    mutation_drivers, amplification_drivers, deletion_drivers : list[str] | str, optional
        Driver genes to profile per alteration type. Default to the curated PDAC driver lists.
    data_lake_path : str, optional
        Directory holding the DepMap files (default: "./data/biomni_data/data_lake").
    mutation_csv_path : str, optional
        Genotype table overriding the default call source (see the synthetic_lethality tools).
    output_csv_path : str, optional
        Write the full landscape table to this CSV path.

    Returns
    -------
    str
        A research log with cohort composition, per-driver group sizes, an alteration frequency
        estimate, a testability verdict per driver and a RECOMMENDED_DRIVERS line for the next stage.

    """
    import pandas as pd

    log = [
        "=" * 78,
        f"PDAC DRIVER LANDSCAPE - testability of each driver in the {cancer_type} cohort",
        f"Run at {datetime.now(tz=_UTC).strftime('%Y-%m-%d %H:%M UTC')}",
        "=" * 78,
        "",
    ]

    try:
        bundle = _load_depmap(data_lake_path)
    except FileNotFoundError as e:
        return f"FAILURE: {e}"

    cohort, match_note = _select_cancer_models(bundle["model"], cancer_type)
    cohort = cohort[cohort["ModelID"].isin(bundle["gene_effect"].index)]
    if len(cohort) == 0:
        return f"FAILURE: no DepMap line with CRISPR data matched cancer_type='{cancer_type}'."

    log.append("STEP 1 | Cohort")
    log.append(f"  {match_note}")
    log.append(f"  Lines with CRISPR gene-effect data: {len(cohort)}")
    subtypes = cohort["OncotreePrimaryDisease"].value_counts().head(5).to_dict()
    log.append(f"  Primary disease composition: {subtypes}")
    log.append("")
    log.append("STEP 2 | Per-driver group sizes (a discovery run needs both arms, not just frequency)")

    rows = []
    plans = [
        ("mutation", _parse_gene_list(mutation_drivers) or PDAC_DRIVERS_MUTATION),
        ("amplification", _parse_gene_list(amplification_drivers) or PDAC_DRIVERS_AMPLIFICATION),
        ("deletion", _parse_gene_list(deletion_drivers) or PDAC_DRIVERS_DELETION),
    ]
    for mode, genes in plans:
        for gene in genes:
            try:
                if mode == "mutation":
                    annotated, source = _annotate_mutation_status(
                        cohort, gene, bundle["data_lake_path"], mutation_csv_path
                    )
                else:
                    annotated, source = _annotate_copy_number_status(
                        cohort, gene, bundle["data_lake_path"], mode, mutation_csv_path
                    )
            except (RuntimeError, ValueError) as e:
                rows.append(
                    {"driver": gene, "alteration": mode, "n_altered": 0, "n_control": 0,
                     "n_unprofiled": len(cohort), "altered_fraction": float("nan"),
                     "verdict": "NO CALLS", "note": str(e)[:120]}
                )
                continue

            n_altered = int((annotated["MutationStatus"] == "MUT").sum())
            n_control = int((annotated["MutationStatus"] == "WT").sum())
            n_unknown = int((annotated["MutationStatus"] == "UNKNOWN").sum())
            profiled = n_altered + n_control
            fraction = n_altered / profiled if profiled else float("nan")

            if min(n_altered, n_control) < MIN_GROUP_SIZE:
                verdict = "UNTESTABLE"
                note = (
                    f"one arm has <{MIN_GROUP_SIZE} lines - a within-PDAC contrast is impossible; "
                    "use a pan-cancer cohort or an isogenic model instead"
                )
            elif min(n_altered, n_control) < INFORMATIVE_GROUP_SIZE:
                verdict = "UNDERPOWERED"
                note = "both arms exist but the smaller is <8 lines; expect no gene to clear a genome-wide FDR"
            else:
                verdict = "TESTABLE"
                note = "both arms are large enough for a genome-wide contrast"
            rows.append(
                {"driver": gene, "alteration": mode, "n_altered": n_altered, "n_control": n_control,
                 "n_unprofiled": n_unknown, "altered_fraction": fraction, "verdict": verdict, "note": note}
            )

    landscape = pd.DataFrame(rows)
    order = {"TESTABLE": 0, "UNDERPOWERED": 1, "UNTESTABLE": 2, "NO CALLS": 3}
    landscape["_order"] = landscape["verdict"].map(order).fillna(4)
    landscape = landscape.sort_values(["_order", "n_altered"], ascending=[True, False]).drop(columns="_order")

    log.append("")
    log.append(f"{'driver':<10}{'alteration':<15}{'n_alt':>7}{'n_ctrl':>8}{'n_NA':>7}{'alt %':>8}  verdict")
    log.append("-" * 78)
    for _, row in landscape.iterrows():
        fraction = "n/a" if pd.isna(row["altered_fraction"]) else f"{row['altered_fraction'] * 100:.0f}%"
        log.append(
            f"{row['driver']:<10}{row['alteration']:<15}{row['n_altered']:>7}{row['n_control']:>8}"
            f"{row['n_unprofiled']:>7}{fraction:>8}  {row['verdict']}"
        )

    testable = landscape[landscape["verdict"] == "TESTABLE"]
    underpowered = landscape[landscape["verdict"] == "UNDERPOWERED"]
    untestable = landscape[landscape["verdict"] == "UNTESTABLE"]

    log.append("")
    log.append("STEP 3 | Reading")
    for _, row in untestable.head(6).iterrows():
        log.append(f"  {row['driver']} ({row['alteration']}): {row['note']}")
    if len(testable) == 0:
        log.append(
            "  No driver gives two adequately sized arms in this cohort. This is the normal situation for "
            "PDAC and it is a property of the cell-line panel, not of the biology: proceed with the "
            "underpowered drivers but treat every output as a ranked hypothesis list, never as an FDR-gated "
            "result, and plan the organoid arm as the real test."
        )

    # What limits a contrast is the SMALLER arm, so recommend on that, and keep the alteration type -
    # SMAD4 mutation and SMAD4 deletion are different experiments on the same gene.
    recommended = pd.concat([testable, underpowered])
    if len(recommended):
        recommended = recommended.assign(
            min_arm=recommended[["n_altered", "n_control"]].min(axis=1)
        ).sort_values("min_arm", ascending=False)
    log.append("")
    log.append(
        "RECOMMENDED_DRIVERS: "
        + ", ".join(
            f"{row['driver']} ({row['alteration']}, smaller arm n={int(row['min_arm'])})"
            for _, row in recommended.iterrows()
        )
    )
    log.append(
        "  (pass one of these as target_mutation to discover_synthetic_lethal_candidates, "
        "discover_allele_resolved_sl_candidates or discover_sl_multichannel)"
    )

    log.append("")
    log.append("QC WARNINGS")
    log.append(
        "  - Alteration fractions here are cell-line frequencies, not patient frequencies. Lines are selected "
        "for growth in 2D culture and under-represent the classical/well-differentiated end of PDAC."
    )
    log.append(
        "  - A TESTABLE verdict means the statistics can run, not that the contrast is biologically clean; "
        "lineage and co-mutation confounding still have to be checked (check_dependency_confounders)."
    )

    if output_csv_path:
        os.makedirs(os.path.dirname(os.path.abspath(output_csv_path)), exist_ok=True)
        landscape.to_csv(output_csv_path, index=False)
        log.append(f"  Full landscape table written to {output_csv_path}")

    log.append("")
    log.append("PROVENANCE")
    for item in bundle["provenance"]:
        log.append(f"  - {item}")
    log.append(f"  - Driver lists: curated TCGA/ICGC PDAC consensus drivers ({len(PDAC_DRIVERS_MUTATION)} mutation, "
               f"{len(PDAC_DRIVERS_AMPLIFICATION)} amplification, {len(PDAC_DRIVERS_DELETION)} deletion)")
    return "\n".join(log)


# ---------------------------------------------------------------------------
# Stage 2: candidate gene -> drug
# ---------------------------------------------------------------------------
def map_sl_candidates_to_drugs(
    candidate_genes,
    data_lake_path: str | None = None,
    min_clinical_phase: str = "Preclinical",
    max_drugs_per_gene: int = 4,
    include_synlethdb_context: bool = True,
    driver_gene: str | None = None,
    output_csv_path: str | None = None,
) -> str:
    """Map synthetic-lethal candidate genes to compounds that actually exist, with clinical phase.

    A candidate that cannot be drugged is not a dead end - it is a CRISPR-only experiment - but the
    distinction has to be made explicitly, before an organoid screen is designed around a target
    with no chemical probe. This tool separates the candidate list into targets with a launched or
    clinical-phase inhibitor, targets with preclinical tool compounds only, and targets with
    nothing, and names the compounds in each case.

    Parameters
    ----------
    candidate_genes : list[str] | str
        Candidate genes, e.g. the CANDIDATE_GENES line from a discovery run.
    data_lake_path : str, optional
        Directory holding broad_repurposing_hub_phase_moa_target_info.parquet.
    min_clinical_phase : str, optional
        Lowest phase to report: "Launched", "Phase 3", "Phase 2", "Phase 1" or "Preclinical"
        (default: "Preclinical", i.e. report everything).
    max_drugs_per_gene : int, optional
        Number of compounds listed per gene, most advanced first (default: 4).
    include_synlethdb_context : bool, optional
        Also report whether the gene is a recorded SynLethDB partner of ``driver_gene`` (default: True).
    driver_gene : str, optional
        Driver used for the SynLethDB context lookup, e.g. "KRAS".
    output_csv_path : str, optional
        Write the full gene-drug table to this CSV path.

    Returns
    -------
    str
        A research log with a per-gene drug table, a druggability verdict per candidate, and
        DRUGGABLE_TARGETS / CRISPR_ONLY_TARGETS lines for the experiment design stage.

    """
    import pandas as pd

    genes = _parse_gene_list(candidate_genes)
    if not genes:
        return "FAILURE: no candidate gene was supplied."

    bundle = _load_repurposing_hub(data_lake_path)
    if bundle is None:
        return (
            "FAILURE: broad_repurposing_hub_phase_moa_target_info.parquet was not found in the data lake. "
            "Drug mapping needs it."
        )

    minimum_rank = CLINICAL_PHASE_RANK.get(min_clinical_phase, 1)
    table = bundle["table"]
    partners = _synlethdb_partners(driver_gene, data_lake_path) if (include_synlethdb_context and driver_gene) else {}

    log = [
        "=" * 78,
        f"DRUG CANDIDATE MAPPING - {len(genes)} synthetic lethal candidates"
        + (f" of {driver_gene.upper()}" if driver_gene else ""),
        f"Run at {datetime.now(tz=_UTC).strftime('%Y-%m-%d %H:%M UTC')}",
        "=" * 78,
        "",
    ]

    rows = []
    clinical, preclinical_only, non_inhibitor_only, undruggable = [], [], [], []
    for gene in genes:
        hits = table[(table["gene"] == gene) & (table["phase_rank"] >= minimum_rank)].copy()
        hits["inhibitory"] = hits["moa"].map(_is_inhibitory)
        # Inhibitory mechanisms first, then clinical phase: a launched substrate is not a target drug.
        hits = hits.sort_values(["inhibitory", "phase_rank"], ascending=[False, False]).drop_duplicates("drug")
        synlethdb_note = ""
        if gene in partners:
            synlethdb_note = f"SynLethDB partner of {driver_gene.upper()} ({', '.join(partners[gene]['sources'])})"

        if len(hits) == 0:
            undruggable.append(gene)
            rows.append({"gene": gene, "drug": "", "clinical_phase": "", "moa": "", "inhibitory": False,
                         "verdict": "CRISPR-ONLY", "synlethdb": synlethdb_note})
            continue

        inhibitors = hits[hits["inhibitory"]]
        if len(inhibitors) == 0:
            verdict = "NO INHIBITOR - annotated compounds are substrates/modulators, not inhibitors"
            non_inhibitor_only.append(gene)
        elif inhibitors.iloc[0]["phase_rank"] >= 2:
            verdict = "CLINICAL-STAGE INHIBITOR"
            clinical.append(gene)
        else:
            verdict = "PRECLINICAL PROBE ONLY"
            preclinical_only.append(gene)
        for _, hit in hits.head(max_drugs_per_gene).iterrows():
            rows.append({"gene": gene, "drug": hit["drug"], "clinical_phase": hit["clinical_phase"],
                         "moa": hit["moa"], "inhibitory": bool(hit["inhibitory"]),
                         "verdict": verdict, "synlethdb": synlethdb_note})

    frame = pd.DataFrame(rows)

    log.append("STEP 1 | Candidate-by-candidate")
    for gene in genes:
        subset = frame[frame["gene"] == gene]
        verdict = subset["verdict"].iloc[0]
        log.append(f"  {gene:<10} {verdict}")
        context = subset["synlethdb"].iloc[0]
        if context:
            log.append(f"             {context}")
        for _, row in subset.iterrows():
            if row["drug"]:
                marker = "" if row["inhibitory"] else "   (non-inhibitory MOA)"
                log.append(f"             - {row['drug']} [{row['clinical_phase']}] {row['moa']}{marker}")
        if verdict == "CRISPR-ONLY":
            log.append("             - no compound in the Repurposing Hub targets this gene")

    log.append("")
    log.append("STEP 2 | Reading for the organoid stage")
    log.append(f"  Clinical-stage inhibitor available : {len(clinical)} ({', '.join(clinical) or 'none'})")
    log.append(f"  Preclinical tool compound only     : {len(preclinical_only)} ({', '.join(preclinical_only) or 'none'})")
    log.append(
        f"  Annotated compounds but no inhibitor: {len(non_inhibitor_only)} "
        f"({', '.join(non_inhibitor_only) or 'none'}) - treat these as CRISPR-only targets"
    )
    log.append(f"  No compound - CRISPR arm required  : {len(undruggable)} ({', '.join(undruggable) or 'none'})")
    log.append(
        "  A drug arm and a CRISPR arm answer different questions: the compound tests whether inhibiting the "
        "protein kills the organoid, the knockout tests whether the gene is required. Discordance between them "
        "is informative (off-target activity, scaffolding function), not a failure of the prediction."
    )

    log.append("")
    log.append(f"DRUGGABLE_TARGETS: {', '.join(clinical + preclinical_only)}")
    log.append(f"CRISPR_ONLY_TARGETS: {', '.join(undruggable + non_inhibitor_only)}")

    log.append("")
    log.append("QC WARNINGS")
    log.append(
        "  - Repurposing Hub target annotations are curated primary targets, not selectivity profiles. "
        "A compound listed for a gene may hit several others at the concentrations an organoid screen uses."
    )
    log.append(
        "  - A target annotation is not a mechanism: cystine is annotated to its transporter SLC7A11 and "
        "glutathione to GPX4. Only compounds with an inhibitory MOA are counted as druggable here, and the "
        "rest are routed to the CRISPR arm."
    )
    log.append(
        "  - Clinical phase describes the compound's development history in its original indication, which "
        "says nothing about activity in PDAC."
    )

    if output_csv_path:
        os.makedirs(os.path.dirname(os.path.abspath(output_csv_path)), exist_ok=True)
        frame.to_csv(output_csv_path, index=False)
        log.append(f"  Gene-drug table written to {output_csv_path} ({len(frame)} rows)")

    log.append("")
    log.append("PROVENANCE")
    log.append(
        f"  - Broad Repurposing Hub: {bundle['path']} ({bundle['n_drugs']} compounds, "
        f"{bundle['n_genes']} annotated target genes)"
    )
    if driver_gene and partners:
        log.append(f"  - SynLethDB context: {len(partners)} recorded partners of {driver_gene.upper()}")
    return "\n".join(log)


# ---------------------------------------------------------------------------
# Stage 3: measured organoid drug response
# ---------------------------------------------------------------------------
def _fit_four_parameter_logistic(log_concentration, viability):
    """Fit viability = lower + (upper - lower) / (1 + 10**((x - log_ic50) * hill)), x = log10(conc).

    ``upper`` is the untreated plateau (low dose), ``lower`` is the maximum-effect plateau (high
    dose, i.e. Emax) and ``log_ic50`` is the curve midpoint between them - the relative EC50, not
    an absolute 50%-of-untreated IC50.

    Returns (parameters, fitted_callable) or (None, None) if the fit does not converge.
    """
    import numpy as np
    from scipy.optimize import curve_fit

    def model(x, upper, lower, log_ic50, hill):
        return lower + (upper - lower) / (1.0 + 10.0 ** ((x - log_ic50) * hill))

    initial = [
        float(np.nanmax(viability)),
        float(np.nanmin(viability)),
        float(np.nanmedian(log_concentration)),
        1.0,
    ]
    bounds = (
        [0.0, -0.5, float(np.nanmin(log_concentration)) - 3, 0.1],
        [2.0, 1.5, float(np.nanmax(log_concentration)) + 3, 10.0],
    )
    try:
        parameters, _ = curve_fit(model, log_concentration, viability, p0=initial, bounds=bounds, maxfev=20000)
    except Exception:
        return None, None
    return parameters, (lambda x: model(x, *parameters))


def analyze_organoid_drug_response(
    screen_csv_path: str,
    organoid_column: str = "organoid",
    drug_column: str = "drug",
    concentration_column: str = "concentration_um",
    viability_column: str = "viability",
    replicate_column: str | None = None,
    genotype_csv_path: str | None = None,
    genotype_column: str = "genotype",
    output_csv_path: str | None = None,
    plot_output_path: str | None = None,
) -> str:
    """Fit measured patient-derived organoid dose-response curves and compare genotype groups.

    Reads the viability readout a PDO drug screen produces, fits a 4-parameter logistic curve per
    (organoid, drug), and reports IC50, AUC and Emax per organoid. When an organoid genotype table
    is supplied, the per-organoid parameters are compared between genotype groups with a Welch
    t-test - that comparison, not the absolute IC50, is what tests a synthetic lethality prediction.

    Expected CSV layout (column names are all parameters, so an existing export can be used as is)::

        organoid,drug,concentration_um,viability,replicate
        PDO-01,MRTX1133,0.01,0.98,1
        PDO-01,MRTX1133,0.1,0.81,1
        ...

    Viability may be a fraction (0-1) or a percentage (0-100); the scale is detected and reported.

    Parameters
    ----------
    screen_csv_path : str
        CSV/TSV of the screen readout.
    organoid_column, drug_column, concentration_column, viability_column : str, optional
        Column names in the screen file.
    replicate_column : str, optional
        Replicate identifier. When given, replicates are kept as separate points in the fit and the
        per-concentration spread is reported.
    genotype_csv_path : str, optional
        CSV mapping organoid to genotype, e.g. "organoid,genotype" with values like "KRAS-G12D" /
        "KRAS-WT". Enables the between-group comparison.
    genotype_column : str, optional
        Genotype column name in that file (default: "genotype").
    output_csv_path : str, optional
        Write the per-organoid, per-drug curve parameters to this CSV.
    plot_output_path : str, optional
        Write a dose-response figure (one panel per drug, curves coloured by genotype) to this PNG.

    Returns
    -------
    str
        A research log with per-curve parameters, genotype-group comparisons, QC on failed or
        censored fits, and a DRUG_RESPONSE_SUMMARY line per drug.

    """
    import numpy as np
    import pandas as pd
    from scipy import stats

    try:
        screen = _read_table(
            screen_csv_path,
            {"organoid": organoid_column, "drug": drug_column,
             "concentration": concentration_column, "viability": viability_column},
            "drug screen",
        )
    except ValueError as e:
        return f"FAILURE: {e}"

    log = [
        "=" * 78,
        "ORGANOID DRUG RESPONSE - measured dose-response, genotype-stratified",
        f"Run at {datetime.now(tz=_UTC).strftime('%Y-%m-%d %H:%M UTC')}",
        "=" * 78,
        "",
    ]

    screen = screen.rename(
        columns={organoid_column: "organoid", drug_column: "drug",
                 concentration_column: "concentration", viability_column: "viability"}
    )
    screen["concentration"] = pd.to_numeric(screen["concentration"], errors="coerce")
    screen["viability"] = pd.to_numeric(screen["viability"], errors="coerce")
    screen = screen.dropna(subset=["concentration", "viability"])
    screen = screen[screen["concentration"] > 0]
    if len(screen) == 0:
        return "FAILURE: no usable row (concentration must be > 0 and viability numeric)."

    scale_note = "fraction (0-1)"
    if screen["viability"].max() > 2:
        screen["viability"] = screen["viability"] / 100.0
        scale_note = "percent (0-100), rescaled to 0-1"

    genotypes = {}
    if genotype_csv_path:
        try:
            metadata = _read_table(
                genotype_csv_path, {"organoid": organoid_column, "genotype": genotype_column}, "genotype table"
            )
            genotypes = dict(zip(metadata[organoid_column], metadata[genotype_column], strict=False))
        except ValueError as e:
            return f"FAILURE: {e}"

    log.append("STEP 1 | Input")
    log.append(f"  Screen file: {screen_csv_path} ({len(screen)} usable measurements)")
    log.append(f"  Organoids: {screen['organoid'].nunique()} | Drugs: {screen['drug'].nunique()}")
    log.append(f"  Viability scale detected: {scale_note}")
    if genotypes:
        counts = pd.Series(list(genotypes.values())).value_counts().to_dict()
        log.append(f"  Genotype groups: {counts}")
    else:
        log.append("  No genotype table supplied - per-organoid parameters only, no group comparison.")

    rows, failures = [], []
    for (organoid, drug), block in screen.groupby(["organoid", "drug"]):
        concentrations = block["concentration"].to_numpy(dtype=float)
        viabilities = block["viability"].to_numpy(dtype=float)
        log_concentration = np.log10(concentrations)
        n_points = len(np.unique(concentrations))
        if n_points < 4:
            failures.append(f"{organoid}/{drug}: only {n_points} distinct concentrations (4 needed for a 4PL fit)")
            continue

        parameters, curve = _fit_four_parameter_logistic(log_concentration, viabilities)
        tested_min, tested_max = float(np.min(concentrations)), float(np.max(concentrations))
        if parameters is None:
            failures.append(f"{organoid}/{drug}: 4PL fit did not converge")
            ic50, emax, censored = np.nan, float(np.min(viabilities)), "fit failed"
            grid = np.linspace(log_concentration.min(), log_concentration.max(), 50)
            means = block.groupby("concentration")["viability"].mean()
            auc = float(np.trapz(means.to_numpy(), np.log10(means.index.to_numpy())) /
                        (np.log10(tested_max) - np.log10(tested_min)))
        else:
            upper, lower, log_ic50, hill = parameters
            ic50 = float(10 ** log_ic50)
            emax = float(lower)  # plateau at the highest dose = maximum achievable effect
            grid = np.linspace(log_concentration.min(), log_concentration.max(), 200)
            auc = float(np.trapz(curve(grid), grid) / (grid[-1] - grid[0]))
            censored = ""
            if ic50 > tested_max:
                censored = f"IC50 > highest tested dose ({tested_max:g})"
            elif ic50 < tested_min:
                censored = f"IC50 < lowest tested dose ({tested_min:g})"

        rows.append(
            {"organoid": organoid, "drug": drug, "genotype": genotypes.get(organoid, ""),
             "ic50_um": ic50, "log10_ic50": float(np.log10(ic50)) if ic50 == ic50 and ic50 > 0 else np.nan,
             "auc": auc, "emax": emax, "n_concentrations": n_points,
             "tested_range_um": f"{tested_min:g}-{tested_max:g}", "censored": censored}
        )

    if not rows:
        log.append("")
        log.append("FAILURE: no curve could be fitted.")
        log.extend(f"  - {item}" for item in failures)
        return "\n".join(log)

    curves = pd.DataFrame(rows)

    log.append("")
    log.append("STEP 2 | Fitted curves (AUC 1.0 = no effect, 0 = complete kill)")
    log.append(f"{'organoid':<12}{'drug':<16}{'genotype':<14}{'IC50 uM':>10}{'AUC':>8}{'Emax':>8}  note")
    log.append("-" * 86)
    for _, row in curves.sort_values(["drug", "auc"]).iterrows():
        ic50_text = "n/a" if not np.isfinite(row["ic50_um"]) else f"{row['ic50_um']:.3g}"
        log.append(
            f"{str(row['organoid']):<12}{str(row['drug']):<16}{str(row['genotype']):<14}"
            f"{ic50_text:>10}{row['auc']:>8.3f}{row['emax']:>8.3f}  {row['censored']}"
        )

    comparisons = []
    if genotypes and curves["genotype"].nunique() >= 2:
        log.append("")
        log.append("STEP 3 | Genotype-group comparison per drug (Welch t-test)")
        for drug, block in curves.groupby("drug"):
            groups = [(name, group) for name, group in block.groupby("genotype") if len(group) >= 2]
            if len(groups) < 2:
                log.append(f"  {drug}: fewer than 2 organoids in at least one genotype group - not tested")
                continue
            (name_a, group_a), (name_b, group_b) = groups[0], groups[1]
            auc_t, auc_p = stats.ttest_ind(group_a["auc"], group_b["auc"], equal_var=False)
            delta_auc = float(group_a["auc"].mean() - group_b["auc"].mean())
            usable_a = group_a["log10_ic50"].dropna()
            usable_b = group_b["log10_ic50"].dropna()
            if len(usable_a) >= 2 and len(usable_b) >= 2:
                _, ic50_p = stats.ttest_ind(usable_a, usable_b, equal_var=False)
                fold = float(10 ** (usable_a.mean() - usable_b.mean()))
                ic50_text = f"IC50 ratio {name_a}/{name_b} = {fold:.2f}x, p={ic50_p:.3f}"
            else:
                ic50_p = np.nan
                ic50_text = "IC50 not comparable (too many censored or failed fits)"
            log.append(
                f"  {drug}: AUC {name_a} {group_a['auc'].mean():.3f} (n={len(group_a)}) vs "
                f"{name_b} {group_b['auc'].mean():.3f} (n={len(group_b)}), delta={delta_auc:+.3f}, p={auc_p:.3f}"
            )
            log.append(f"      {ic50_text}")
            comparisons.append(
                {"drug": drug, "group_a": name_a, "group_b": name_b, "mean_auc_a": float(group_a["auc"].mean()),
                 "mean_auc_b": float(group_b["auc"].mean()), "delta_auc": delta_auc, "auc_p": float(auc_p),
                 "ic50_p": float(ic50_p) if ic50_p == ic50_p else np.nan,
                 "n_a": len(group_a), "n_b": len(group_b)}
            )

        for record in comparisons:
            direction = "more sensitive" if record["delta_auc"] < 0 else "less sensitive"
            log.append(
                f"DRUG_RESPONSE_SUMMARY: {record['drug']} | {record['group_a']} is {direction} than "
                f"{record['group_b']} | delta_AUC={record['delta_auc']:+.3f} | p={record['auc_p']:.3f} | "
                f"n={record['n_a']}v{record['n_b']}"
            )

    log.append("")
    log.append("QC WARNINGS")
    for item in failures[:10]:
        log.append(f"  - {item}")
    censored_rows = curves[curves["censored"] != ""]
    if len(censored_rows):
        log.append(
            f"  - {len(censored_rows)} curve(s) have an IC50 outside the tested concentration range. Those IC50 "
            "values are extrapolations; the AUC over the tested range is the honest summary and is what the "
            "group comparison uses."
        )
    small_groups = curves.groupby("genotype").size()
    if genotypes and (small_groups < 3).any():
        log.append(
            f"  - Genotype groups are small {small_groups.to_dict()}; a p-value from 2-3 organoids per arm is "
            "descriptive. Organoid-to-organoid variability in PDAC is large."
        )
    log.append(
        "  - Viability at a single timepoint cannot separate cytotoxicity from cytostasis; a drug that holds "
        "AUC at 0.5 by arresting growth and one that kills half the cells look identical here."
    )

    if plot_output_path:
        try:
            _plot_dose_response(screen, curves, genotypes, plot_output_path)
            log.append(f"  Dose-response figure written to {plot_output_path}")
        except Exception as e:
            log.append(f"  Dose-response figure FAILED: {e}")

    if output_csv_path:
        os.makedirs(os.path.dirname(os.path.abspath(output_csv_path)), exist_ok=True)
        curves.to_csv(output_csv_path, index=False)
        log.append(f"  Curve parameter table written to {output_csv_path} ({len(curves)} curves)")

    log.append("")
    log.append("PROVENANCE")
    log.append(f"  - Screen data: {screen_csv_path} ({len(screen)} measurements, {scale_note})")
    if genotype_csv_path:
        log.append(f"  - Genotypes: {genotype_csv_path}")
    log.append(
        "  - Model: 4-parameter logistic on log10 concentration (scipy curve_fit, bounded). IC50 is the "
        "fitted curve midpoint (relative EC50), Emax the fitted plateau at the highest dose, and AUC the mean "
        "fitted viability over the tested log-concentration range"
    )
    return "\n".join(log)


def _plot_dose_response(screen, curves, genotypes: dict, output_path: str) -> None:
    """One panel per drug: measured points and fitted curves, coloured by genotype."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np
    import seaborn as sns

    sns.set_theme(style="whitegrid", context="notebook")
    drugs = sorted(screen["drug"].unique())
    n_columns = min(3, len(drugs))
    n_rows = int(np.ceil(len(drugs) / n_columns))
    figure, axes = plt.subplots(n_rows, n_columns, figsize=(4.6 * n_columns, 4.0 * n_rows), squeeze=False)

    group_names = sorted({str(g) for g in genotypes.values()}) if genotypes else []
    palette = dict(zip(group_names, sns.color_palette("Set1", max(len(group_names), 1)), strict=False))

    for index, drug in enumerate(drugs):
        axis = axes[index // n_columns][index % n_columns]
        block = screen[screen["drug"] == drug]
        for organoid, organoid_block in block.groupby("organoid"):
            genotype = str(genotypes.get(organoid, ""))
            color = palette.get(genotype, "0.4")
            means = organoid_block.groupby("concentration")["viability"].mean()
            axis.scatter(means.index, means.to_numpy(), s=22, color=color, alpha=0.8)
            row = curves[(curves["organoid"] == organoid) & (curves["drug"] == drug)]
            if len(row) and np.isfinite(row.iloc[0]["ic50_um"]):
                grid = np.linspace(np.log10(means.index.min()), np.log10(means.index.max()), 100)
                parameters, curve = _fit_four_parameter_logistic(
                    np.log10(organoid_block["concentration"].to_numpy(dtype=float)),
                    organoid_block["viability"].to_numpy(dtype=float),
                )
                if curve is not None:
                    axis.plot(10 ** grid, curve(grid), color=color, linewidth=1.6,
                              label=f"{organoid} ({genotype})" if genotype else str(organoid))
        axis.set_xscale("log")
        axis.axhline(0.5, color="0.6", linestyle="--", linewidth=0.9)
        axis.set_title(drug, fontsize=11)
        axis.set_xlabel("concentration (uM, log scale)")
        axis.set_ylabel("viability (fraction of control)" if index % n_columns == 0 else "")
        axis.legend(fontsize=7, loc="lower left", frameon=True)

    for empty in range(len(drugs), n_rows * n_columns):
        axes[empty // n_columns][empty % n_columns].axis("off")
    figure.suptitle("Patient-derived organoid dose-response", fontsize=13, y=0.995)
    figure.tight_layout(rect=(0, 0, 1, 0.98))
    os.makedirs(os.path.dirname(os.path.abspath(output_path)), exist_ok=True)
    figure.savefig(output_path, dpi=200, bbox_inches="tight")
    plt.close(figure)


# ---------------------------------------------------------------------------
# Stage 4: measured organoid CRISPR validation
# ---------------------------------------------------------------------------
def analyze_crispr_validation(
    validation_csv_path: str,
    organoid_column: str = "organoid",
    gene_column: str = "gene",
    viability_column: str = "viability",
    control_label: str = "NTC",
    replicate_column: str | None = "replicate",
    genotype_csv_path: str | None = None,
    genotype_column: str = "genotype",
    mutant_genotype_label: str | None = None,
    output_csv_path: str | None = None,
) -> str:
    """Turn measured organoid CRISPR knockout viability into per-gene effects and a selectivity test.

    Each knockout is normalised to the non-targeting control of the same organoid, so organoid-level
    differences in growth rate cannot masquerade as a knockout effect. Two questions are then
    answered separately, because a prediction of synthetic lethality requires both:

    1. does the knockout reduce viability at all (per organoid, vs its own control), and
    2. is the reduction larger in driver-mutant organoids than in wild-type ones (the selectivity
       that distinguishes a synthetic lethal partner from a core fitness gene).

    Expected CSV layout::

        organoid,gene,viability,replicate
        PDO-01,NTC,1.02,1
        PDO-01,TEAD1,0.61,1
        ...

    Parameters
    ----------
    validation_csv_path : str
        CSV/TSV of the validation readout.
    organoid_column, gene_column, viability_column : str, optional
        Column names in the validation file.
    control_label : str, optional
        Value in the gene column marking the non-targeting control (default: "NTC").
    replicate_column : str, optional
        Replicate identifier; replicates are what the per-knockout t-test is computed over.
    genotype_csv_path : str, optional
        CSV mapping organoid to genotype. Required for the selectivity test.
    genotype_column : str, optional
        Genotype column name in that file (default: "genotype").
    mutant_genotype_label : str, optional
        Which genotype value counts as driver-mutant. Defaults to the most frequent value that is
        not obviously wild type (any label containing "WT" or "wild" is treated as the control).
    output_csv_path : str, optional
        Write the per-gene, per-organoid effect table to this CSV.

    Returns
    -------
    str
        A research log with per-organoid knockout effects, a per-gene selectivity verdict and
        CRISPR_VALIDATION lines that ``compare_prediction_with_experiment`` can consume.

    """
    import numpy as np
    import pandas as pd
    from scipy import stats

    try:
        data = _read_table(
            validation_csv_path,
            {"organoid": organoid_column, "gene": gene_column, "viability": viability_column},
            "CRISPR validation",
        )
    except ValueError as e:
        return f"FAILURE: {e}"

    data = data.rename(columns={organoid_column: "organoid", gene_column: "gene", viability_column: "viability"})
    data["viability"] = pd.to_numeric(data["viability"], errors="coerce")
    data = data.dropna(subset=["viability"])
    if data["viability"].max() > 2:
        data["viability"] = data["viability"] / 100.0

    log = [
        "=" * 78,
        "ORGANOID CRISPR VALIDATION - knockout effect and genotype selectivity",
        f"Run at {datetime.now(tz=_UTC).strftime('%Y-%m-%d %H:%M UTC')}",
        "=" * 78,
        "",
    ]

    controls = data[data["gene"].astype(str).str.upper() == control_label.upper()]
    if len(controls) == 0:
        return (
            f"FAILURE: no row with gene == '{control_label}'. The non-targeting control defines the baseline; "
            "pass the right `control_label`."
        )

    genotypes = {}
    if genotype_csv_path:
        try:
            metadata = _read_table(
                genotype_csv_path, {"organoid": organoid_column, "genotype": genotype_column}, "genotype table"
            )
            genotypes = dict(zip(metadata[organoid_column], metadata[genotype_column], strict=False))
        except ValueError as e:
            return f"FAILURE: {e}"

    log.append("STEP 1 | Input")
    log.append(f"  Validation file: {validation_csv_path} ({len(data)} measurements)")
    log.append(f"  Organoids: {data['organoid'].nunique()} | Knockouts: {data['gene'].nunique() - 1} + control")
    log.append(f"  Control label: {control_label} ({len(controls)} measurements)")

    rows = []
    for (organoid, gene), block in data.groupby(["organoid", "gene"]):
        if str(gene).upper() == control_label.upper():
            continue
        control_block = controls[controls["organoid"] == organoid]["viability"]
        if len(control_block) == 0:
            continue
        control_mean = float(control_block.mean())
        if control_mean <= 0:
            continue
        relative = block["viability"].to_numpy(dtype=float) / control_mean
        mean_relative = float(np.mean(relative))
        log2_fold_change = float(np.log2(max(mean_relative, 1e-6)))
        if len(block) >= 2 and len(control_block) >= 2:
            _, p_value = stats.ttest_ind(block["viability"], control_block, equal_var=False)
        else:
            p_value = np.nan
        rows.append(
            {"organoid": organoid, "gene": gene, "genotype": genotypes.get(organoid, ""),
             "relative_viability": mean_relative, "log2_fold_change": log2_fold_change,
             "p_vs_control": float(p_value) if p_value == p_value else np.nan,
             "n_replicates": len(block), "n_control_replicates": len(control_block)}
        )

    if not rows:
        return "FAILURE: no knockout could be normalised to a control in the same organoid."
    effects = pd.DataFrame(rows)

    log.append("")
    log.append("STEP 2 | Knockout effect per organoid (relative to that organoid's own control)")
    log.append(f"{'organoid':<12}{'gene':<12}{'genotype':<14}{'rel. viab':>10}{'log2FC':>9}{'p':>10}{'n':>4}")
    log.append("-" * 72)
    for _, row in effects.sort_values(["gene", "log2_fold_change"]).iterrows():
        p_text = "n/a" if not np.isfinite(row["p_vs_control"]) else f"{row['p_vs_control']:.3f}"
        log.append(
            f"{str(row['organoid']):<12}{str(row['gene']):<12}{str(row['genotype']):<14}"
            f"{row['relative_viability']:>10.3f}{row['log2_fold_change']:>9.3f}{p_text:>10}{row['n_replicates']:>4}"
        )

    summary_rows = []
    if genotypes and effects["genotype"].nunique() >= 2:
        labels = sorted(effects["genotype"].dropna().unique())
        control_labels = [label for label in labels if "WT" in str(label).upper() or "WILD" in str(label).upper()]
        if mutant_genotype_label:
            mutant_label = mutant_genotype_label
            wildtype_label = next((label for label in labels if label != mutant_label), None)
        elif control_labels:
            wildtype_label = control_labels[0]
            mutant_label = next((label for label in labels if label != wildtype_label), None)
        else:
            mutant_label, wildtype_label = labels[0], labels[1]

        log.append("")
        log.append(f"STEP 3 | Genotype selectivity: {mutant_label} vs {wildtype_label} organoids")
        log.append(
            f"{'gene':<12}{'log2FC mut':>12}{'log2FC wt':>11}{'delta':>9}{'p':>9}{'n':>8}  verdict"
        )
        log.append("-" * 76)
        for gene, block in effects.groupby("gene"):
            mutant_block = block[block["genotype"] == mutant_label]["log2_fold_change"]
            wildtype_block = block[block["genotype"] == wildtype_label]["log2_fold_change"]
            if len(mutant_block) == 0 or len(wildtype_block) == 0:
                continue
            delta = float(mutant_block.mean() - wildtype_block.mean())
            if len(mutant_block) >= 2 and len(wildtype_block) >= 2:
                _, p_value = stats.ttest_ind(mutant_block, wildtype_block, equal_var=False)
            else:
                p_value = np.nan
            depleted = mutant_block.mean() <= -0.5
            if depleted and delta <= -0.5 and (not np.isfinite(p_value) or p_value < 0.05):
                verdict = "SELECTIVE - supports synthetic lethality"
            elif depleted and delta > -0.5:
                verdict = "NON-SELECTIVE - kills both genotypes (core fitness)"
            elif not depleted:
                verdict = "NO EFFECT - knockout does not reduce viability"
            else:
                verdict = "TREND ONLY - direction right, not significant"
            p_text = "n/a" if not np.isfinite(p_value) else f"{p_value:.3f}"
            log.append(
                f"{str(gene):<12}{mutant_block.mean():>12.3f}{wildtype_block.mean():>11.3f}{delta:>9.3f}"
                f"{p_text:>9}{f'{len(mutant_block)}v{len(wildtype_block)}':>8}  {verdict}"
            )
            summary_rows.append(
                {"gene": gene, "log2fc_mutant": float(mutant_block.mean()),
                 "log2fc_wildtype": float(wildtype_block.mean()), "selectivity_delta": delta,
                 "selectivity_p": float(p_value) if p_value == p_value else np.nan,
                 "n_mutant": len(mutant_block), "n_wildtype": len(wildtype_block), "verdict": verdict}
            )

        log.append("")
        for record in summary_rows:
            log.append(
                f"CRISPR_VALIDATION: {record['gene']} | delta={record['selectivity_delta']:+.3f} | "
                f"p={record['selectivity_p'] if record['selectivity_p'] == record['selectivity_p'] else float('nan'):.3f} | "
                f"{record['verdict']}"
            )
    else:
        log.append("")
        log.append(
            "STEP 3 | Genotype selectivity not tested: supply `genotype_csv_path` with at least two genotype "
            "groups. Without it, a knockout that kills every organoid is indistinguishable from a synthetic "
            "lethal one."
        )

    log.append("")
    log.append("QC WARNINGS")
    thin = effects[effects["n_replicates"] < 3]
    if len(thin):
        log.append(f"  - {len(thin)} knockout/organoid pairs have <3 replicates; their p-values are unstable.")
    log.append(
        "  - Viability normalised to a non-targeting control does not control for cutting toxicity; an "
        "amplified locus gives a growth defect from DNA damage alone. A multi-sgRNA design is what separates "
        "that from a real gene effect."
    )
    log.append(
        "  - A selective knockout effect in organoids is evidence for the dependency, not for a therapeutic "
        "window: the matched normal-organoid arm is what tests toxicity."
    )

    if output_csv_path:
        os.makedirs(os.path.dirname(os.path.abspath(output_csv_path)), exist_ok=True)
        table = pd.DataFrame(summary_rows) if summary_rows else effects
        table.to_csv(output_csv_path, index=False)
        log.append(f"  Table written to {output_csv_path} ({len(table)} rows)")

    log.append("")
    log.append("PROVENANCE")
    log.append(f"  - Validation data: {validation_csv_path}")
    if genotype_csv_path:
        log.append(f"  - Genotypes: {genotype_csv_path}")
    log.append("  - Statistics: Welch t-test on replicates (knockout vs control) and on log2FC between genotypes")
    return "\n".join(log)


# ---------------------------------------------------------------------------
# Stage 5: did the prediction hold?
# ---------------------------------------------------------------------------
def compare_prediction_with_experiment(
    prediction_csv_path: str,
    experiment_csv_path: str,
    prediction_gene_column: str = "gene",
    prediction_score_column: str = "effect_difference",
    experiment_gene_column: str = "gene",
    experiment_effect_column: str = "selectivity_delta",
    prediction_lower_is_stronger: bool = True,
    experiment_lower_is_stronger: bool = True,
    experiment_hit_threshold: float = -0.5,
    top_k: int = 5,
    n_permutations: int = 10000,
    output_csv_path: str | None = None,
) -> str:
    """Score Biomni's predictions against what the organoid experiment measured.

    This is the step that makes the project falsifiable. It joins the predicted ranking to the
    measured effects, then reports: rank concordance over the tested genes, precision@k with a
    permutation p-value (so a hit rate is read against what random selection of the same number of
    genes would give), and - printed first - the predictions the experiment falsified.

    A small tested set is the normal case here, and the tool refuses to dress it up: with fewer than
    5 shared genes no correlation is computed, only the per-gene confirm/falsify table.

    Parameters
    ----------
    prediction_csv_path : str
        CSV written by a discovery/ranking tool (e.g. output_csv_path of
        stratify_ovarian_cancer_dependency_by_mutation or discover_synthetic_lethal_candidates).
    experiment_csv_path : str
        CSV written by ``analyze_crispr_validation`` (or any table of measured per-gene effects;
        for a drug screen, map drugs back to their target gene first).
    prediction_gene_column, prediction_score_column : str, optional
        Columns in the prediction table (default: "gene", "effect_difference").
    experiment_gene_column, experiment_effect_column : str, optional
        Columns in the experiment table (default: "gene", "selectivity_delta").
    prediction_lower_is_stronger, experiment_lower_is_stronger : bool, optional
        Direction conventions. Both default to True (more negative = stronger dependency).
    experiment_hit_threshold : float, optional
        Measured effect at or beyond which a gene counts as an experimental hit (default: -0.5,
        i.e. a halving of viability in log2 units).
    top_k : int, optional
        Size of the predicted top set used for precision@k (default: 5).
    n_permutations : int, optional
        Permutations for the precision@k null (default: 10000).
    output_csv_path : str, optional
        Write the joined prediction/experiment table to this CSV.

    Returns
    -------
    str
        A research log with the falsified predictions first, then confirmed ones, rank concordance,
        precision@k against a permutation null, and an overall VERDICT line.

    """
    import numpy as np
    import pandas as pd
    from scipy import stats

    try:
        predictions = _read_table(
            prediction_csv_path,
            {"gene": prediction_gene_column, "score": prediction_score_column}, "prediction table",
        )
        experiments = _read_table(
            experiment_csv_path,
            {"gene": experiment_gene_column, "effect": experiment_effect_column}, "experiment table",
        )
    except ValueError as e:
        return f"FAILURE: {e}"

    log = [
        "=" * 78,
        "PREDICTION vs EXPERIMENT - did the predicted vulnerabilities hold in organoids?",
        f"Run at {datetime.now(tz=_UTC).strftime('%Y-%m-%d %H:%M UTC')}",
        "=" * 78,
        "",
    ]

    predictions = predictions[[prediction_gene_column, prediction_score_column]].copy()
    predictions.columns = ["gene", "predicted_score"]
    predictions["gene"] = predictions["gene"].astype(str).str.upper().str.strip()
    predictions = predictions.dropna(subset=["predicted_score"]).drop_duplicates("gene")
    # Rank 1 = most strongly predicted, whichever direction the score runs in.
    predictions["predicted_rank"] = predictions["predicted_score"].rank(
        ascending=prediction_lower_is_stronger, method="min"
    ).astype(int)

    experiments = experiments[[experiment_gene_column, experiment_effect_column]].copy()
    experiments.columns = ["gene", "measured_effect"]
    experiments["gene"] = experiments["gene"].astype(str).str.upper().str.strip()
    experiments = experiments.dropna(subset=["measured_effect"]).drop_duplicates("gene")

    joined = predictions.merge(experiments, on="gene", how="inner")
    if len(joined) == 0:
        return (
            "FAILURE: no gene is present in both tables. Check that the gene symbols match "
            f"(prediction e.g. {predictions['gene'].head(3).tolist()}, "
            f"experiment e.g. {experiments['gene'].head(3).tolist()})."
        )

    if experiment_lower_is_stronger:
        joined["experimental_hit"] = joined["measured_effect"] <= experiment_hit_threshold
    else:
        joined["experimental_hit"] = joined["measured_effect"] >= experiment_hit_threshold
    joined = joined.sort_values("predicted_rank")

    tested_k = min(top_k, len(joined))
    predicted_top = joined.head(tested_k)
    n_hits_total = int(joined["experimental_hit"].sum())
    n_hits_in_top = int(predicted_top["experimental_hit"].sum())
    precision_at_k = n_hits_in_top / tested_k if tested_k else float("nan")
    base_rate = n_hits_total / len(joined)

    log.append("STEP 1 | Overlap")
    log.append(f"  Predicted genes: {len(predictions)} | Experimentally measured genes: {len(experiments)}")
    log.append(f"  Genes in both (the only ones this comparison can speak about): {len(joined)}")
    log.append(
        f"  Experimental hits (measured effect {'<=' if experiment_lower_is_stronger else '>='} "
        f"{experiment_hit_threshold}): {n_hits_total} of {len(joined)}"
    )

    falsified = predicted_top[~predicted_top["experimental_hit"]]
    confirmed = predicted_top[predicted_top["experimental_hit"]]
    missed = joined[(joined["experimental_hit"]) & (joined["predicted_rank"] > tested_k)]

    log.append("")
    log.append(f"STEP 2 | FALSIFIED predictions (predicted top {tested_k}, not confirmed experimentally)")
    if len(falsified) == 0:
        log.append("  none - every gene in the predicted top set reached the experimental hit threshold")
    for _, row in falsified.iterrows():
        log.append(
            f"  {row['gene']:<10} predicted rank {int(row['predicted_rank']):<3} "
            f"(score {row['predicted_score']:+.3f}) -> measured {row['measured_effect']:+.3f} "
            f"(threshold {experiment_hit_threshold})"
        )

    log.append("")
    log.append("STEP 3 | CONFIRMED predictions")
    if len(confirmed) == 0:
        log.append("  none")
    for _, row in confirmed.iterrows():
        log.append(
            f"  {row['gene']:<10} predicted rank {int(row['predicted_rank']):<3} "
            f"(score {row['predicted_score']:+.3f}) -> measured {row['measured_effect']:+.3f}"
        )

    log.append("")
    log.append("STEP 4 | MISSED by the prediction (experimental hits ranked outside the top set)")
    if len(missed) == 0:
        log.append("  none")
    for _, row in missed.iterrows():
        log.append(
            f"  {row['gene']:<10} predicted rank {int(row['predicted_rank']):<3} -> "
            f"measured {row['measured_effect']:+.3f}"
        )

    log.append("")
    log.append("STEP 5 | Quantitative concordance")
    log.append(f"  precision@{tested_k} = {n_hits_in_top}/{tested_k} = {precision_at_k:.2f}")
    log.append(f"  base rate (hits among all tested genes) = {n_hits_total}/{len(joined)} = {base_rate:.2f}")

    permutation_p = float("nan")
    if tested_k >= len(joined):
        log.append(
            "  permutation p not computable: the predicted top set IS the entire tested set, so there is no "
            "null to compare it against. A falsifiable test has to carry some genes the method ranked LOW into "
            "the experiment as negative controls, or use a smaller top_k."
        )
    elif n_hits_total == 0:
        log.append("  permutation p not computable: no gene reached the experimental hit threshold.")
    if n_hits_total > 0 and len(joined) > tested_k:
        rng = np.random.default_rng(0)
        hit_flags = joined["experimental_hit"].to_numpy()
        null_counts = np.array(
            [rng.permutation(hit_flags)[:tested_k].sum() for _ in range(n_permutations)]
        )
        permutation_p = float((null_counts >= n_hits_in_top).mean())
        log.append(
            f"  permutation p (random {tested_k} genes from the tested set reach {n_hits_in_top} hits) "
            f"= {permutation_p:.4f}"
        )

    if len(joined) >= 5:
        spearman = stats.spearmanr(joined["predicted_score"], joined["measured_effect"])
        sign = 1 if prediction_lower_is_stronger == experiment_lower_is_stronger else -1
        log.append(
            f"  Spearman rho (predicted score vs measured effect) = {spearman.statistic:.3f}, "
            f"p = {spearman.pvalue:.4f}"
        )
        log.append(
            f"  A concordant method gives rho {'> 0' if sign > 0 else '< 0'} here, given the direction "
            "conventions supplied."
        )
    else:
        log.append(
            f"  Rank correlation not computed: only {len(joined)} shared genes. With n < 5 a correlation "
            "coefficient is noise, and the per-gene table above is the whole result."
        )

    log.append("")
    enrichment_testable = np.isfinite(permutation_p)
    if precision_at_k >= 0.5 and enrichment_testable and permutation_p < 0.05:
        verdict = "PREDICTIVE - the predicted top set is enriched for experimental hits beyond chance"
    elif precision_at_k >= 0.5 and not enrichment_testable:
        verdict = (
            "CONSISTENT BUT UNTESTED - most of the predicted top set was confirmed, but with no low-ranked "
            "genes in the experiment the enrichment cannot be separated from the base rate"
        )
    elif n_hits_in_top > 0:
        verdict = "WEAK - some predictions confirmed, but not beyond what picking genes at random would give"
    else:
        verdict = "NOT PREDICTIVE - no gene in the predicted top set reached the experimental threshold"
    log.append(f"VERDICT: {verdict}")
    permutation_text = f"{permutation_p:.4f}" if enrichment_testable else "not computable"
    log.append(
        f"OVERALL_CONCORDANCE: precision@{tested_k}={precision_at_k:.2f} | base_rate={base_rate:.2f} | "
        f"permutation_p={permutation_text} | n_shared={len(joined)}"
    )

    log.append("")
    log.append("QC WARNINGS")
    log.append(
        "  - Only genes carried into the experiment can be scored. Candidates dropped before the wet lab are "
        "neither confirmed nor falsified, and a precision computed over a set the prediction itself chose is "
        "optimistic: it never sees the genes the method ranked low."
    )
    if len(joined) < 10:
        log.append(
            f"  - {len(joined)} shared genes is a small test set; precision@k moves in steps of "
            f"{1 / tested_k:.2f} and one organoid failure changes the verdict."
        )
    log.append(
        "  - The hit threshold is a pre-specified convention, not a biological boundary. Report it with the "
        "result and do not tune it after seeing the data."
    )

    if output_csv_path:
        os.makedirs(os.path.dirname(os.path.abspath(output_csv_path)), exist_ok=True)
        joined.to_csv(output_csv_path, index=False)
        log.append(f"  Joined table written to {output_csv_path} ({len(joined)} genes)")

    log.append("")
    log.append("PROVENANCE")
    log.append(f"  - Predictions: {prediction_csv_path} (column {prediction_score_column})")
    log.append(f"  - Experiment: {experiment_csv_path} (column {experiment_effect_column})")
    log.append(
        f"  - Hit threshold {experiment_hit_threshold}; precision@{tested_k} null from {n_permutations} "
        "permutations of the measured hit labels"
    )
    return "\n".join(log)


# ---------------------------------------------------------------------------
# Stage 6: HTML report
# ---------------------------------------------------------------------------
def _literature_evidence(
    disease: str,
    driver: str,
    genes: list,
    max_papers: int = 6,
    email: str | None = None,
    api_key: str | None = None,
) -> dict:
    """Fetch real PubMed evidence per candidate gene. Never invents a citation.

    Returns ``{gene: {"score", "interpretation", "n_papers", "query", "papers": [...]}}``; a gene
    with no PubMed hit gets an empty paper list rather than a plausible-looking one.
    """
    import time

    from biomni.tool.synthetic_lethality import _entrez_efetch, _entrez_esearch, _score_literature

    evidence = {}
    delay = 0.4 if api_key else 0.75
    for gene in genes:
        pair = f"({driver}[Title/Abstract] AND {gene}[Title/Abstract])"
        contextual = f"{pair} AND ({disease}[Title/Abstract] OR pancreatic[Title/Abstract])"
        queries = [("disease-context", contextual), ("gene-pair", pair)]

        pmids: list[str] = []
        used = []
        for label, query in queries:
            hits = _entrez_esearch(query, max_papers, email, api_key)
            time.sleep(delay)
            new = [pmid for pmid in hits if pmid not in pmids]
            used.append(f"[{label}] {len(hits)} hit(s)")
            pmids.extend(new)
            if len(pmids) >= max_papers:
                break

        records = _entrez_efetch(pmids[:max_papers], email, api_key)
        time.sleep(delay)
        scored = _score_literature(records, gene, driver, disease)

        # A PubMed hit for "KRAS AND <gene> AND pancreatic" is often a panel paper that merely lists
        # the gene. Count the records whose text really carries BOTH symbols, so a high volume score
        # built on co-mention cannot pass for evidence about the interaction.
        def mentions_pair(record):
            text = f"{record.get('title', '')} {record.get('abstract', '')}".upper()
            return driver.upper() in text and gene.upper() in text

        paired = [record for record in records if mentions_pair(record)]
        evidence[gene] = {
            "score": scored["score"],
            "interpretation": scored["interpretation"],
            "n_papers": scored["n_papers"],
            "n_pair_papers": len(paired),
            "refuting_pmids": scored.get("refuting_pmids", []),
            "query": "; ".join(used),
            "papers": [
                {
                    "pmid": record.get("pmid", ""),
                    "title": record.get("title", ""),
                    "journal": record.get("journal", ""),
                    "year": record.get("year"),
                    "pair": mentions_pair(record),
                }
                for record in records
                if record.get("pmid")
            ],
        }
    return evidence


def _read_optional_csv(path: str):
    import pandas as pd

    if path and os.path.exists(path):
        return pd.read_csv(path)
    return None


def generate_pdac_report(
    run_dir: str = "./pdac_run",
    output_html_path: str | None = None,
    disease: str = "Pancreatic cancer",
    driver: str = "KRAS",
    include_literature: bool = True,
    max_papers_per_gene: int = 6,
    top_candidates: int = 10,
    standalone: bool = True,
    email: str | None = None,
    api_key: str | None = None,
) -> str:
    """Build a self-contained HTML report from the outputs of the PDAC pipeline.

    Reads whatever the run produced - driver landscape, DepMap predictions, drug mapping, and (when
    the wet-lab stages have run) organoid curves, CRISPR validation and the concordance table - and
    renders them as charts plus the tables behind them. References are fetched live from PubMed per
    candidate gene, so every citation in the report is a real PMID with a link; genes with no
    literature are shown as having none rather than being given a plausible citation.

    Parameters
    ----------
    run_dir : str, optional
        Directory holding the pipeline CSVs (default: "./pdac_run").
    output_html_path : str, optional
        Where to write the HTML (default: ``<run_dir>/pdac_report.html``).
    disease, driver : str, optional
        Context used for the PubMed queries and the report header.
    include_literature : bool, optional
        Query PubMed for per-candidate references (default: True; needs internet).
    max_papers_per_gene : int, optional
        PubMed records retrieved per candidate (default: 6).
    top_candidates : int, optional
        Number of candidates charted in detail (default: 10).
    standalone : bool, optional
        True writes a complete HTML document for local viewing; False writes a body fragment
        (title + style + content) for embedding in a host that supplies the document skeleton.
    email, api_key : str, optional
        NCBI Entrez contact e-mail and API key.

    Returns
    -------
    str
        A short log naming the report path, the sections rendered and the evidence counts.

    """
    import json

    landscape = _read_optional_csv(os.path.join(run_dir, "driver_landscape.csv"))
    predictions = _read_optional_csv(os.path.join(run_dir, "prediction_depmap.csv"))
    drugs = _read_optional_csv(os.path.join(run_dir, "drug_candidates.csv"))
    curves = _read_optional_csv(os.path.join(run_dir, "pdo_curves.csv"))
    crispr = _read_optional_csv(os.path.join(run_dir, "crispr_summary.csv"))
    concordance = _read_optional_csv(os.path.join(run_dir, "concordance.csv"))

    if landscape is None and predictions is None:
        return (
            f"FAILURE: neither driver_landscape.csv nor prediction_depmap.csv was found in {run_dir}. "
            "Run `python run_pdac_agent.py --mode direct --stages 1-6` first."
        )

    output_html_path = output_html_path or os.path.join(run_dir, "pdac_report.html")
    payload: dict = {
        "generated_at": datetime.now(tz=_UTC).strftime("%Y-%m-%d %H:%M UTC"),
        "disease": disease,
        "driver": driver,
        "run_dir": os.path.abspath(run_dir),
        "sections": [],
    }

    if landscape is not None:
        landscape = landscape.copy()
        landscape["min_arm"] = landscape[["n_altered", "n_control"]].min(axis=1)
        payload["landscape"] = landscape.sort_values("min_arm", ascending=False).to_dict("records")
        payload["sections"].append("driver_landscape")

    candidates: list = []
    if predictions is not None:
        charted = predictions.head(top_candidates).copy()
        candidates = charted["gene"].astype(str).str.upper().tolist()
        payload["predictions"] = charted.to_dict("records")
        payload["n_predictions"] = int(len(predictions))
        payload["sections"].append("predictions")

    if drugs is not None:
        payload["drugs"] = drugs.to_dict("records")
        payload["sections"].append("drugs")

    for name, frame in (("curves", curves), ("crispr", crispr), ("concordance", concordance)):
        if frame is not None:
            payload[name] = frame.to_dict("records")
            payload["sections"].append(name)

    literature = {}
    if include_literature and candidates:
        literature = _literature_evidence(
            disease, driver, candidates, max_papers=max_papers_per_gene, email=email, api_key=api_key
        )
        payload["sections"].append("literature")
    payload["literature"] = literature

    document = _REPORT_TEMPLATE.replace("__REPORT_DATA__", json.dumps(payload, default=str))
    if standalone:
        document = (
            "<!doctype html>\n<html lang=\"ko\">\n<head>\n<meta charset=\"utf-8\">\n"
            "<meta name=\"viewport\" content=\"width=device-width, initial-scale=1\">\n"
            + document
            + "\n</html>\n"
        )

    os.makedirs(os.path.dirname(os.path.abspath(output_html_path)), exist_ok=True)
    with open(output_html_path, "w", encoding="utf-8") as handle:
        handle.write(document)

    n_papers = sum(len(record["papers"]) for record in literature.values())
    log = [
        f"HTML report written to {output_html_path} ({os.path.getsize(output_html_path) / 1024:.0f} KB)",
        f"  Sections rendered: {', '.join(payload['sections']) or 'none'}",
        f"  Drivers profiled: {len(payload.get('landscape', []))}",
        f"  Candidates charted: {len(candidates)} of {payload.get('n_predictions', 0)}",
        f"  References: {n_papers} PubMed records across {len(literature)} candidates "
        f"({sum(1 for r in literature.values() if not r['papers'])} candidate(s) with no PubMed hit)",
    ]
    if concordance is None:
        log.append(
            "  Experimental stages not present: the report shows the prediction side only and states the "
            "pre-specified criteria the organoid experiment has to meet."
        )
    return "\n".join(log)


# The report shell. Charts are inline SVG drawn from the embedded payload: no external script or
# stylesheet is loaded except the font, so the file works offline and inside a strict CSP.
_REPORT_TEMPLATE = r"""<title>PDAC Vulnerability Dossier</title>
<link rel="preconnect" href="https://fonts.googleapis.com">
<link rel="preconnect" href="https://fonts.gstatic.com" crossorigin>
<link rel="stylesheet" href="https://fonts.googleapis.com/css2?family=IBM+Plex+Mono:wght@400;500&family=IBM+Plex+Sans+Condensed:wght@500;600;700&family=IBM+Plex+Serif:ital,wght@0,400;0,600;1,400&display=swap">
<style>
:root{
  color-scheme: light;
  --page:#f4f6f5; --surface:#fcfdfc; --surface-2:#eef1f0;
  --ink:#0e1211; --ink-2:#4c5553; --muted:#828b88;
  --rule:#dde2e0; --grid:#e6eae8; --axis:#c3c9c6;
  --accent:#1c5cab; --series-1:#2a78d6; --series-2:#eb6834;
  --good:#0ca30c; --warn:#fab219; --crit:#d03b3b;
  --shadow:0 1px 2px rgba(14,18,17,.05), 0 8px 24px -16px rgba(14,18,17,.25);
  --sans:"IBM Plex Sans Condensed", system-ui, -apple-system, "Segoe UI", sans-serif;
  --serif:"IBM Plex Serif", Georgia, "Times New Roman", serif;
  --mono:"IBM Plex Mono", ui-monospace, SFMono-Regular, Menlo, monospace;
}
@media (prefers-color-scheme: dark){
  :root:where(:not([data-theme="light"])){
    color-scheme: dark;
    --page:#0e100f; --surface:#1a1a19; --surface-2:#222423;
    --ink:#ffffff; --ink-2:#c3c2b7; --muted:#898781;
    --rule:#2c2e2c; --grid:#2c2c2a; --axis:#383835;
    --accent:#6da7ec; --series-1:#3987e5; --series-2:#d95926;
    --good:#0ca30c; --warn:#fab219; --crit:#d03b3b;
    --shadow:0 1px 2px rgba(0,0,0,.4), 0 8px 24px -16px rgba(0,0,0,.8);
  }
}
:root[data-theme="dark"]{
  color-scheme: dark;
  --page:#0e100f; --surface:#1a1a19; --surface-2:#222423;
  --ink:#ffffff; --ink-2:#c3c2b7; --muted:#898781;
  --rule:#2c2e2c; --grid:#2c2c2a; --axis:#383835;
  --accent:#6da7ec; --series-1:#3987e5; --series-2:#d95926;
  --good:#0ca30c; --warn:#fab219; --crit:#d03b3b;
  --shadow:0 1px 2px rgba(0,0,0,.4), 0 8px 24px -16px rgba(0,0,0,.8);
}
*{box-sizing:border-box}
body{margin:0;background:var(--page);color:var(--ink);font-family:var(--serif);
  font-size:16px;line-height:1.62;-webkit-font-smoothing:antialiased}
a{color:var(--accent);text-underline-offset:2px}
a:focus-visible,summary:focus-visible{outline:2px solid var(--accent);outline-offset:2px;border-radius:2px}
.wrap{max-width:1080px;margin:0 auto;padding:0 28px 96px}
header.masthead{border-bottom:1px solid var(--rule);padding:56px 0 28px;margin-bottom:8px}
.eyebrow{font-family:var(--sans);font-size:12px;font-weight:600;letter-spacing:.14em;
  text-transform:uppercase;color:var(--muted)}
h1{font-family:var(--sans);font-weight:700;font-size:clamp(30px,4.4vw,46px);line-height:1.08;
  letter-spacing:-.015em;margin:14px 0 10px;text-wrap:balance;max-width:20ch}
.question{font-family:var(--serif);font-style:italic;color:var(--ink-2);max-width:62ch;margin:0 0 20px}
.runmeta{font-family:var(--mono);font-size:12px;color:var(--muted);display:flex;flex-wrap:wrap;gap:6px 18px}
.tiles{display:grid;grid-template-columns:repeat(auto-fit,minmax(170px,1fr));gap:12px;margin:28px 0 12px}
.tile{background:var(--surface);border:1px solid var(--rule);border-radius:3px;padding:14px 16px;box-shadow:var(--shadow)}
.tile .k{font-family:var(--sans);font-size:11px;font-weight:600;letter-spacing:.1em;text-transform:uppercase;color:var(--muted)}
.tile .v{font-family:var(--sans);font-size:30px;font-weight:700;line-height:1.1;margin-top:6px;letter-spacing:-.02em}
.tile .n{font-family:var(--mono);font-size:11.5px;color:var(--ink-2);margin-top:4px;line-height:1.45}
main{display:grid;grid-template-columns:72px minmax(0,1fr);gap:0 28px}
section.stage{grid-column:1 / -1;display:grid;grid-template-columns:subgrid;padding:40px 0 8px;border-top:1px solid var(--rule)}
section.stage:first-of-type{border-top:none}
.rail{font-family:var(--mono);font-size:12px;color:var(--muted);padding-top:6px}
.rail .num{display:block;font-family:var(--sans);font-size:26px;font-weight:700;color:var(--axis);line-height:1}
.rail .state{display:inline-block;margin-top:8px;font-size:10.5px;letter-spacing:.06em;text-transform:uppercase}
.rail .state.done{color:var(--good)} .rail .state.pending{color:var(--muted)}
.body{min-width:0}
h2{font-family:var(--sans);font-weight:600;font-size:23px;letter-spacing:-.01em;margin:0 0 6px;text-wrap:balance}
.finding{margin:0 0 22px;max-width:66ch;color:var(--ink-2)}
.finding strong{color:var(--ink);font-weight:600}
figure{margin:0 0 18px;background:var(--surface);border:1px solid var(--rule);border-radius:3px;
  padding:18px 18px 12px;box-shadow:var(--shadow)}
figcaption{font-family:var(--mono);font-size:11.5px;color:var(--muted);margin-top:10px;line-height:1.5}
.legend{display:flex;flex-wrap:wrap;gap:8px 16px;font-family:var(--sans);font-size:12px;
  font-weight:500;color:var(--ink-2);margin-bottom:12px}
.legend i{display:inline-block;width:10px;height:10px;border-radius:2px;margin-right:6px;vertical-align:-1px}
svg{display:block;width:100%;height:auto;overflow:visible}
svg text{font-family:var(--mono);font-size:11px;fill:var(--ink-2)}
svg text.lab{font-family:var(--sans);font-weight:600;font-size:12px;fill:var(--ink)}
svg text.val{font-family:var(--mono);font-size:11px;fill:var(--ink-2)}
svg .gridline{stroke:var(--grid);stroke-width:1}
svg .axis{stroke:var(--axis);stroke-width:1}
svg .hit{fill:transparent;cursor:crosshair}
details{margin:0 0 18px;border:1px solid var(--rule);border-radius:3px;background:var(--surface)}
summary{font-family:var(--sans);font-size:12.5px;font-weight:600;letter-spacing:.04em;
  text-transform:uppercase;color:var(--ink-2);padding:11px 14px;cursor:pointer}
.tablewrap{overflow-x:auto;padding:0 14px 14px}
table{border-collapse:collapse;width:100%;font-family:var(--mono);font-size:12px;
  font-variant-numeric:tabular-nums;white-space:nowrap}
th{text-align:left;font-family:var(--sans);font-size:11px;letter-spacing:.06em;text-transform:uppercase;
  color:var(--muted);font-weight:600;padding:6px 12px 6px 0;border-bottom:1px solid var(--rule)}
td{padding:5px 12px 5px 0;border-bottom:1px solid var(--grid);color:var(--ink-2)}
td.gene{color:var(--ink);font-weight:500}
.chip{display:inline-flex;align-items:center;gap:5px;font-family:var(--sans);font-size:11px;
  font-weight:600;letter-spacing:.03em;padding:1px 7px;border-radius:2px;border:1px solid currentColor}
.chip.good{color:var(--good)} .chip.warn{color:#9a6b00} .chip.crit{color:var(--crit)}
:root[data-theme="dark"] .chip.warn,
:root:where(:not([data-theme="light"])) .chip.warn{color:var(--warn)}
@media (prefers-color-scheme: light){:root:not([data-theme="dark"]) .chip.warn{color:#9a6b00}}
.refs{list-style:none;padding:0;margin:0}
.refs li{padding:11px 0;border-bottom:1px solid var(--grid);max-width:74ch}
.refs .gene{font-family:var(--sans);font-weight:700;font-size:13px;letter-spacing:.04em}
.refs .meta{font-family:var(--mono);font-size:11.5px;color:var(--muted)}
.refs .title{display:block;margin-top:2px}
.none{font-family:var(--mono);font-size:12px;color:var(--muted)}
.callout{border-left:3px solid var(--accent);background:var(--surface-2);padding:14px 16px;
  margin:0 0 18px;border-radius:0 3px 3px 0}
.callout p{margin:0;max-width:64ch;font-size:15px}
.callout .h{font-family:var(--sans);font-size:11.5px;font-weight:700;letter-spacing:.1em;
  text-transform:uppercase;color:var(--accent);margin-bottom:5px}
footer{border-top:1px solid var(--rule);margin-top:44px;padding-top:20px;font-family:var(--mono);
  font-size:11.5px;color:var(--muted);line-height:1.7}
#tip{position:fixed;z-index:20;pointer-events:none;opacity:0;transition:opacity .12s;
  background:var(--ink);color:var(--page);font-family:var(--mono);font-size:11.5px;line-height:1.5;
  padding:7px 9px;border-radius:3px;max-width:280px;box-shadow:var(--shadow)}
@media (prefers-reduced-motion: reduce){*{transition:none!important;animation:none!important}}
@media (max-width:720px){
  main{grid-template-columns:1fr}
  section.stage{display:block}
  .rail{display:flex;gap:12px;align-items:baseline;padding-bottom:10px}
  .rail .num{font-size:20px}
}
</style>

<div id="tip" role="status" aria-live="polite"></div>
<div class="wrap">
  <header class="masthead">
    <div class="eyebrow">Biomni &middot; translational pipeline run</div>
    <h1>PDAC Vulnerability Dossier</h1>
    <p class="question">췌장암의 유전적 취약성과 합성치사 관계를 발굴하고, 예측된 약물 효과를 환자유래 오가노이드에서 기능적으로 검증할 수 있는가?</p>
    <div class="runmeta" id="runmeta"></div>
    <div class="tiles" id="tiles"></div>
  </header>
  <main id="stages"></main>
  <footer id="foot"></footer>
</div>

<script>
const REPORT = __REPORT_DATA__;

/* ---------- small helpers ---------- */
const NS = "http://www.w3.org/2000/svg";
const tip = document.getElementById("tip");
const fmt = (v, d = 2) => (v === null || v === undefined || Number.isNaN(Number(v))) ? "n/a" : Number(v).toFixed(d);
const sci = v => (v === null || v === undefined || v === "" || Number.isNaN(Number(v))) ? "n/a"
  : (Math.abs(Number(v)) < 0.001 ? Number(v).toExponential(1) : Number(v).toFixed(3));
const esc = s => String(s ?? "").replace(/[&<>"]/g, c => ({ "&": "&amp;", "<": "&lt;", ">": "&gt;", '"': "&quot;" }[c]));
function svgEl(tag, attrs) {
  const node = document.createElementNS(NS, tag);
  for (const [key, value] of Object.entries(attrs || {})) node.setAttribute(key, value);
  return node;
}
function hover(node, html) {
  node.addEventListener("pointermove", event => {
    tip.innerHTML = html;
    tip.style.opacity = 1;
    const box = tip.getBoundingClientRect();
    let x = event.clientX + 14, y = event.clientY + 14;
    if (x + box.width > window.innerWidth - 8) x = event.clientX - box.width - 14;
    if (y + box.height > window.innerHeight - 8) y = event.clientY - box.height - 14;
    tip.style.left = x + "px"; tip.style.top = y + "px";
  });
  node.addEventListener("pointerleave", () => { tip.style.opacity = 0; });
}

/* ---------- chart: horizontal bars ---------- */
function barChart(rows, options) {
  const opt = Object.assign({ width: 760, rowHeight: 26, labelWidth: 108, valueWidth: 118, max: null }, options);
  const plotWidth = opt.width - opt.labelWidth - opt.valueWidth;
  const height = rows.length * opt.rowHeight + 26;
  const svg = svgEl("svg", { viewBox: `0 0 ${opt.width} ${height}`, role: "img" });
  const max = opt.max || Math.max(...rows.map(r => r.value), 1);
  const ticks = 4;
  for (let i = 0; i <= ticks; i++) {
    const x = opt.labelWidth + (plotWidth * i) / ticks;
    svg.appendChild(svgEl("line", { class: "gridline", x1: x, y1: 0, x2: x, y2: rows.length * opt.rowHeight }));
    const label = svgEl("text", { x: x, y: height - 8, "text-anchor": i === 0 ? "start" : "middle" });
    label.textContent = opt.tickFormat ? opt.tickFormat(max * i / ticks) : (max * i / ticks).toFixed(opt.tickDigits ?? 0);
    svg.appendChild(label);
  }
  svg.appendChild(svgEl("line", { class: "axis", x1: opt.labelWidth, y1: 0, x2: opt.labelWidth, y2: rows.length * opt.rowHeight }));
  rows.forEach((row, index) => {
    const y = index * opt.rowHeight;
    const barHeight = 13;
    const width = Math.max(2, (row.value / max) * plotWidth);
    const name = svgEl("text", { class: "lab", x: opt.labelWidth - 10, y: y + barHeight + 1, "text-anchor": "end" });
    name.textContent = row.label;
    svg.appendChild(name);
    svg.appendChild(svgEl("rect", {
      x: opt.labelWidth, y: y + 4, width: width, height: barHeight, rx: 3, fill: row.color || "var(--series-1)"
    }));
    const value = svgEl("text", { class: "val", x: opt.labelWidth + width + 8, y: y + barHeight + 1 });
    value.textContent = row.valueLabel;
    svg.appendChild(value);
    const hit = svgEl("rect", { class: "hit", x: 0, y: y, width: opt.width, height: opt.rowHeight });
    hover(hit, row.tip);
    svg.appendChild(hit);
  });
  return svg;
}

/* ---------- chart: paired dots (selectivity) ---------- */
function dumbbell(rows, options) {
  const opt = Object.assign({ width: 760, rowHeight: 26, labelWidth: 108 }, options);
  const plotWidth = opt.width - opt.labelWidth - 46;
  const height = rows.length * opt.rowHeight + 26;
  const svg = svgEl("svg", { viewBox: `0 0 ${opt.width} ${height}`, role: "img" });
  for (let i = 0; i <= 4; i++) {
    const x = opt.labelWidth + (plotWidth * i) / 4;
    svg.appendChild(svgEl("line", { class: "gridline", x1: x, y1: 0, x2: x, y2: rows.length * opt.rowHeight }));
    const label = svgEl("text", { x: x, y: height - 8, "text-anchor": i === 0 ? "start" : "middle" });
    label.textContent = (i * 25) + "%";
    svg.appendChild(label);
  }
  svg.appendChild(svgEl("line", { class: "axis", x1: opt.labelWidth, y1: 0, x2: opt.labelWidth, y2: rows.length * opt.rowHeight }));
  rows.forEach((row, index) => {
    const y = index * opt.rowHeight + 11;
    const xa = opt.labelWidth + (row.a / 100) * plotWidth;
    const xb = opt.labelWidth + (row.b / 100) * plotWidth;
    const name = svgEl("text", { class: "lab", x: opt.labelWidth - 10, y: y + 4, "text-anchor": "end" });
    name.textContent = row.label;
    svg.appendChild(name);
    svg.appendChild(svgEl("line", { x1: xa, y1: y, x2: xb, y2: y, stroke: "var(--axis)", "stroke-width": 2 }));
    svg.appendChild(svgEl("circle", { cx: xb, cy: y, r: 5, fill: "var(--series-2)", stroke: "var(--surface)", "stroke-width": 2 }));
    svg.appendChild(svgEl("circle", { cx: xa, cy: y, r: 5, fill: "var(--series-1)", stroke: "var(--surface)", "stroke-width": 2 }));
    const hit = svgEl("rect", { class: "hit", x: 0, y: index * opt.rowHeight, width: opt.width, height: opt.rowHeight });
    hover(hit, row.tip);
    svg.appendChild(hit);
  });
  return svg;
}

/* ---------- chart: segmented bar ---------- */
function segmentBar(segments, total) {
  const width = 760, height = 54;
  const svg = svgEl("svg", { viewBox: `0 0 ${width} ${height}`, role: "img" });
  let x = 0;
  segments.filter(s => s.value > 0).forEach(segment => {
    const w = (segment.value / total) * width;
    svg.appendChild(svgEl("rect", { x: x, y: 6, width: Math.max(0, w - 2), height: 26, rx: 3, fill: segment.color }));
    if (w > 34) {
      const label = svgEl("text", { class: "val", x: x + 7, y: 24, fill: "#fcfdfc" });
      label.textContent = segment.value;
      svg.appendChild(label);
    }
    const hit = svgEl("rect", { class: "hit", x: x, y: 0, width: w, height: height });
    hover(hit, `<b>${esc(segment.label)}</b><br>${segment.value} / ${total} candidates<br>${esc(segment.genes || "")}`);
    svg.appendChild(hit);
    x += w;
  });
  return svg;
}

/* ---------- page assembly ---------- */
function tile(key, value, note) {
  return `<div class="tile"><div class="k">${esc(key)}</div><div class="v">${esc(value)}</div><div class="n">${note}</div></div>`;
}
function table(columns, rows) {
  const head = columns.map(c => `<th>${esc(c.title)}</th>`).join("");
  const body = rows.map(row => "<tr>" + columns.map(c => `<td class="${c.cls || ""}">${c.render(row)}</td>`).join("") + "</tr>").join("");
  return `<details><summary>Data table</summary><div class="tablewrap"><table><thead><tr>${head}</tr></thead><tbody>${body}</tbody></table></div></details>`;
}
function stage(number, state, title, finding, nodes) {
  const section = document.createElement("section");
  section.className = "stage";
  section.innerHTML = `<div class="rail"><span class="num">${number}</span>
    <span class="state ${state === "done" ? "done" : "pending"}">${state === "done" ? "● 완료" : "○ 대기"}</span></div>
    <div class="body"><h2>${esc(title)}</h2><div class="finding">${finding}</div></div>`;
  const body = section.querySelector(".body");
  nodes.forEach(node => body.appendChild(node));
  document.getElementById("stages").appendChild(section);
  return body;
}
function figure(svg, caption, legend) {
  const wrap = document.createElement("figure");
  if (legend) {
    const box = document.createElement("div");
    box.className = "legend";
    box.innerHTML = legend;
    wrap.appendChild(box);
  }
  wrap.appendChild(svg);
  const cap = document.createElement("figcaption");
  cap.innerHTML = caption;
  wrap.appendChild(cap);
  return wrap;
}
function html(markup) {
  const node = document.createElement("div");
  node.innerHTML = markup;
  return node;
}

document.getElementById("runmeta").innerHTML =
  [`생성 ${esc(REPORT.generated_at)}`, `드라이버 ${esc(REPORT.driver)}`, `맥락 ${esc(REPORT.disease)}`,
   `출력 ${esc(REPORT.run_dir)}`].map(s => `<span>${s}</span>`).join("");

const landscape = REPORT.landscape || [];
const predictions = REPORT.predictions || [];
const drugs = REPORT.drugs || [];
const literature = REPORT.literature || {};
const crispr = REPORT.crispr || [];
const concordance = REPORT.concordance || [];
const curves = REPORT.curves || [];

const driverRow = landscape.find(r => String(r.driver).toUpperCase() === String(REPORT.driver).toUpperCase() && r.alteration === "mutation");
const testable = landscape.filter(r => r.verdict === "TESTABLE");
const topGene = predictions.length ? predictions[0].gene : "n/a";
const litScores = Object.values(literature).map(v => v.score);
const nRefs = Object.values(literature).reduce((sum, v) => sum + v.papers.length, 0);

document.getElementById("tiles").innerHTML = [
  tile("검정 가능한 드라이버", `${testable.length} / ${landscape.length}`, "양 군 모두 통계가 가능한 조합"),
  driverRow ? tile(`${REPORT.driver} 변이 세포주`, `${driverRow.n_altered} vs ${driverRow.n_control}`,
    `야생형 ${driverRow.n_control}주 — ${esc(driverRow.verdict)}`) : "",
  tile("합성치사 후보", String(REPORT.n_predictions || predictions.length), `FDR 통과, 최상위 ${esc(topGene)}`),
  tile("PubMed 참고문헌", String(nRefs), `${Object.keys(literature).length}개 후보에 대한 실제 PMID`)
].join("");

/* Stage 1 - driver landscape */
if (landscape.length) {
  const statusColor = v => v === "TESTABLE" ? "var(--good)" : v === "UNDERPOWERED" ? "var(--warn)" : "var(--crit)";
  const rows = landscape.slice(0, 14).map(r => ({
    label: `${r.driver} ${r.alteration.slice(0, 3)}`,
    value: r.min_arm,
    valueLabel: `${r.n_altered}v${r.n_control}`,
    color: statusColor(r.verdict),
    tip: `<b>${esc(r.driver)}</b> (${esc(r.alteration)})<br>변이군 ${r.n_altered} · 대조군 ${r.n_control} · 미프로파일 ${r.n_unprofiled}<br>빈도 ${(r.altered_fraction * 100).toFixed(0)}%<br><b>${esc(r.verdict)}</b>`
  }));
  const body = stage("01", "done", "어떤 드라이버가 실제로 검정 가능한가",
    `세포주 코호트에서 대비를 결정하는 것은 <strong>작은 쪽 군의 크기</strong>입니다. ${REPORT.driver}는 빈도가 가장 높지만(${driverRow ? (driverRow.altered_fraction * 100).toFixed(0) : "?"}%) 바로 그 이유로 야생형 대조군이 ${driverRow ? driverRow.n_control : "?"}주뿐이라 <strong>UNDERPOWERED</strong>입니다. 검정력이 확보되는 조합은 ${testable.map(r => r.driver + " " + r.alteration).slice(0, 4).join(", ")}입니다.`,
    [figure(barChart(rows, { tickDigits: 0 }),
      "막대 길이 = 두 군 중 작은 쪽의 세포주 수(검정력을 결정하는 값). 우측 숫자는 변이군 vs 대조군.",
      `<span><i style="background:var(--good)"></i>TESTABLE (≥8)</span><span><i style="background:var(--warn)"></i>UNDERPOWERED (3–7)</span><span><i style="background:var(--crit)"></i>UNTESTABLE (&lt;3)</span>`)]);
  body.appendChild(html(table([
    { title: "Driver", cls: "gene", render: r => esc(r.driver) },
    { title: "Alteration", render: r => esc(r.alteration) },
    { title: "n altered", render: r => r.n_altered },
    { title: "n control", render: r => r.n_control },
    { title: "Freq", render: r => (r.altered_fraction * 100).toFixed(0) + "%" },
    { title: "Verdict", render: r => `<span class="chip ${r.verdict === "TESTABLE" ? "good" : r.verdict === "UNDERPOWERED" ? "warn" : "crit"}">${esc(r.verdict)}</span>` }
  ], landscape)));
}

/* Stage 2 - predictions */
if (predictions.length) {
  const effectRows = predictions.map(p => ({
    label: p.gene,
    value: Math.abs(p.effect_difference),
    valueLabel: Number(p.effect_difference).toFixed(2),
    color: "var(--series-1)",
    tip: `<b>${esc(p.gene)}</b><br>효과 차이 ${fmt(p.effect_difference, 3)} (변이 ${fmt(p.mutant_mean_effect, 2)} vs 야생형 ${fmt(p.wildtype_mean_effect, 2)})<br>p=${sci(p.p_value)} · q=${sci(p.q_value)}<br>변이주 의존 ${fmt(p.pct_mutant_dependent, 0)}% · 전체 ${fmt(p.pct_all_lines_dependent, 0)}%`
  }));
  const selRows = predictions.map(p => ({
    label: p.gene,
    a: Number(p.pct_mutant_dependent),
    b: Number(p.pct_all_lines_dependent),
    tip: `<b>${esc(p.gene)}</b><br>변이주 의존 ${fmt(p.pct_mutant_dependent, 0)}%<br>전체 세포주 의존 ${fmt(p.pct_all_lines_dependent, 0)}%<br>차이가 클수록 선택적 취약점, 작으면 공통 필수 유전자`
  }));
  const body = stage("02", "done", "합성치사 후보와 그 선택성",
    `DepMap CRISPR 스크린에서 ${REPORT.driver} 변이군이 야생형보다 유의하게 더 의존하는 유전자입니다. 양성 대조가 통과했습니다 — ${REPORT.driver} 자신이 자기 변이주에서 가장 깊게 떨어집니다. <strong>효과 크기만으로는 부족합니다</strong>: 전체 세포주에서도 똑같이 필수인 유전자는 합성치사가 아니라 공통 필수 유전자이므로, 아래 두 번째 차트의 두 점이 벌어진 후보만 오가노이드에 투입할 가치가 있습니다.`,
    [
      figure(barChart(effectRows, { tickDigits: 2 }),
        "효과 차이 = 변이군 평균 − 야생형 평균 gene effect. 음수일수록 변이군이 더 큰 타격을 받습니다(막대는 절댓값).", ""),
      figure(dumbbell(selRows),
        "두 점의 간격이 선택성입니다. 겹쳐 있으면 변이와 무관하게 필수인 유전자입니다.",
        `<span><i style="background:var(--series-1)"></i>${esc(REPORT.driver)} 변이주 의존 비율</span><span><i style="background:var(--series-2)"></i>전체 세포주 의존 비율</span>`)
    ]);
  body.appendChild(html(table([
    { title: "Rank", render: p => p.rank },
    { title: "Gene", cls: "gene", render: p => esc(p.gene) },
    { title: "MUT mean", render: p => fmt(p.mutant_mean_effect, 3) },
    { title: "WT mean", render: p => fmt(p.wildtype_mean_effect, 3) },
    { title: "Δ", render: p => fmt(p.effect_difference, 3) },
    { title: "p", render: p => sci(p.p_value) },
    { title: "q", render: p => sci(p.q_value) },
    { title: "%MUT dep", render: p => fmt(p.pct_mutant_dependent, 0) },
    { title: "%all dep", render: p => fmt(p.pct_all_lines_dependent, 0) }
  ], predictions)));
}

/* Stage 3 - druggability */
if (drugs.length) {
  const byGene = new Map();
  drugs.forEach(row => {
    if (!byGene.has(row.gene)) byGene.set(row.gene, { verdict: row.verdict, drugs: [] });
    if (row.drug) byGene.get(row.gene).drugs.push(row);
  });
  const classify = verdict => String(verdict).startsWith("CLINICAL") ? "clinical"
    : String(verdict).startsWith("PRECLINICAL") ? "preclinical"
    : String(verdict).startsWith("NO INHIBITOR") ? "noninhibitor" : "crispr";
  const buckets = { clinical: [], preclinical: [], noninhibitor: [], crispr: [] };
  byGene.forEach((value, gene) => buckets[classify(value.verdict)].push(gene));
  const segments = [
    { label: "임상단계 저해제 있음", value: buckets.clinical.length, color: "var(--series-1)", genes: buckets.clinical.join(", ") },
    { label: "전임상 probe만", value: buckets.preclinical.length, color: "var(--series-2)", genes: buckets.preclinical.join(", ") },
    { label: "화합물은 있으나 저해제 아님", value: buckets.noninhibitor.length, color: "var(--warn)", genes: buckets.noninhibitor.join(", ") },
    { label: "화합물 없음 — CRISPR 전용", value: buckets.crispr.length, color: "var(--axis)", genes: buckets.crispr.join(", ") }
  ];
  const total = byGene.size;
  const body = stage("03", "done", "약물로 공략 가능한 후보와 CRISPR 전용 후보",
    `Broad Repurposing Hub 매핑입니다. 표적 주석(annotation)은 기전이 아니므로 <strong>MOA가 저해(inhibitor/antagonist/degrader)인 화합물만</strong> 약물 arm으로 셉니다 — 그러지 않으면 기질인 cystine이 SLC7A11의 "저해제"로 잡힙니다. 화합물이 없는 후보는 실패가 아니라 <strong>실험 설계의 분기</strong>입니다: 약물 arm 대신 CRISPR arm으로 검증합니다.`,
    [figure(segmentBar(segments, total),
      "후보 " + total + "개의 공략 가능성 분류. 각 구간에 마우스를 올리면 해당 유전자 목록이 보입니다.",
      segments.map(s => `<span><i style="background:${s.color}"></i>${esc(s.label)} (${s.value})</span>`).join(""))]);
  const drugRows = drugs.filter(d => d.drug);
  body.appendChild(html(table([
    { title: "Gene", cls: "gene", render: d => esc(d.gene) },
    { title: "Drug", render: d => esc(d.drug) },
    { title: "Phase", render: d => esc(d.clinical_phase) },
    { title: "MOA", render: d => esc(d.moa) },
    { title: "Inhibitory", render: d => (d.inhibitory === true || d.inhibitory === "True") ? "yes" : "no" }
  ], drugRows)));
}

/* Stage 4 - literature */
if (Object.keys(literature).length) {
  const genes = Object.keys(literature);
  const rows = genes.map(gene => ({
    label: gene,
    value: literature[gene].score,
    valueLabel: `${literature[gene].score}  (쌍 ${literature[gene].n_pair_papers ?? 0}/${literature[gene].n_papers})`,
    color: literature[gene].score >= 65 ? "var(--good)" : literature[gene].score >= 35 ? "var(--series-1)" : "var(--axis)",
    tip: `<b>${esc(gene)}</b><br>${esc(literature[gene].interpretation)}<br>논문 ${literature[gene].n_papers}편 · 점수 ${literature[gene].score}/100<br>${esc(REPORT.driver)}와 함께 다룬 논문 ${literature[gene].n_pair_papers ?? 0}편`
  })).sort((a, b) => b.value - a.value);
  const noHit = genes.filter(g => literature[g].papers.length === 0);
  const body = stage("04", "done", "문헌 근거",
    `각 후보에 대해 PubMed를 <strong>실제로 조회한</strong> 결과입니다(NCBI Entrez). 점수는 논문 수, 유전자 쌍 동시출현, 합성치사 관련 표현, 직접 실험 근거, 질환 맥락, 최신성에서 반박 표현 패널티를 뺀 값입니다. ${noHit.length ? `<strong>${noHit.join(", ")}</strong>는 PubMed 근거가 없습니다 — 없는 것을 없다고 표시하며, 그럴듯한 인용으로 채우지 않습니다.` : ""} 문헌이 적다는 것은 틀렸다는 뜻이 아니라 <strong>아직 검증되지 않았다</strong>는 뜻이므로, 통계 근거와 독립적으로 읽어야 합니다.`,
    [figure(barChart(rows, { max: 100, tickDigits: 0 }),
      "0–100 문헌 <em>기록량</em> 점수이지 진실의 양이 아닙니다. 막대 옆 괄호는 검색된 논문 중 드라이버와 후보를 <em>함께</em> 다룬 편수 — 이 값이 0이면 높은 점수라도 단순 동시 언급입니다(예: 유전자 패널 논문).", "")]);
  const list = document.createElement("ul");
  list.className = "refs";
  genes.forEach(gene => {
    const entry = literature[gene];
    if (!entry.papers.length) {
      list.appendChild(html(`<li><span class="gene">${esc(gene)}</span> <span class="none">— PubMed 검색 결과 없음 (${esc(entry.query)})</span></li>`).firstChild);
      return;
    }
    entry.papers.forEach(paper => {
      list.appendChild(html(`<li><span class="gene">${esc(gene)}</span>
        <span class="meta">PMID ${esc(paper.pmid)} · ${esc(paper.journal || "")} ${paper.year || ""}</span>
        <span class="title">${paper.pair ? "" : '<span class="none">[단순 언급] </span>'}${esc(paper.title)} <a href="https://pubmed.ncbi.nlm.nih.gov/${esc(paper.pmid)}/" target="_blank" rel="noopener">PubMed ↗</a></span></li>`).firstChild);
    });
  });
  const details = document.createElement("details");
  details.innerHTML = `<summary>참고문헌 ${nRefs}건 (PMID 링크)</summary>`;
  const pad = document.createElement("div");
  pad.className = "tablewrap";
  pad.appendChild(list);
  details.appendChild(pad);
  body.appendChild(details);
}

/* Stage 5 - experimental validation */
{
  const hasExperiment = crispr.length || curves.length;
  const nodes = [];
  if (curves.length) {
    const rows = curves.map(c => ({
      label: `${c.organoid} ${String(c.drug).slice(0, 10)}`,
      value: Number(c.auc),
      valueLabel: Number(c.auc).toFixed(2),
      color: "var(--series-1)",
      tip: `<b>${esc(c.organoid)}</b> · ${esc(c.drug)}<br>유전형 ${esc(c.genotype || "n/a")}<br>IC50 ${fmt(c.ic50_um, 3)} µM · AUC ${fmt(c.auc, 3)} · Emax ${fmt(c.emax, 3)}`
    }));
    nodes.push(figure(barChart(rows, { max: 1, tickDigits: 1, labelWidth: 150 }),
      "AUC 1.0 = 효과 없음, 0 = 완전 사멸. 오가노이드별 용량반응 곡선에서 계산.", ""));
  }
  if (crispr.length) {
    const rows = crispr.map(c => ({
      label: c.gene,
      value: Math.abs(Number(c.selectivity_delta)),
      valueLabel: Number(c.selectivity_delta).toFixed(2),
      color: String(c.verdict).startsWith("SELECTIVE") ? "var(--good)" : String(c.verdict).startsWith("NON-SELECTIVE") ? "var(--warn)" : "var(--axis)",
      tip: `<b>${esc(c.gene)}</b><br>변이 오가노이드 log2FC ${fmt(c.log2fc_mutant, 2)} · 야생형 ${fmt(c.log2fc_wildtype, 2)}<br>선택성 Δ ${fmt(c.selectivity_delta, 3)} · p=${sci(c.selectivity_p)}<br>${esc(c.verdict)}`
    }));
    nodes.push(figure(barChart(rows, { tickDigits: 2 }),
      "선택성 Δ = 변이 오가노이드 log2FC − 야생형 log2FC. 초록은 선택적, 주황은 두 유전형 모두 사멸(공통 필수).", ""));
  }
  if (!hasExperiment) {
    nodes.push(html(`<div class="callout"><div class="h">사전 지정 기준</div>
      <p>오가노이드 데이터가 아직 투입되지 않았습니다. 예측을 반증 가능하게 만들려면, 측정 전에 다음을 고정해야 합니다 —
      ① 적중 기준: 변이 오가노이드에서 log2FC ≤ −0.5 <em>그리고</em> 야생형 대비 선택성 Δ ≤ −0.5,
      ② KO 패널에 <strong>예측이 낮게 평가한 유전자도 포함</strong>(없으면 precision@k가 base rate와 같아져 우연과 구분 불가),
      ③ 정상 췌장 오가노이드 대조군으로 치료 창(window) 확인.</p></div>`));
  }
  stage("05", crispr.length || curves.length ? "done" : "pending", "오가노이드 기능 검증",
    crispr.length || curves.length
      ? `측정된 오가노이드 데이터입니다. 약물 arm은 유전형군 간 AUC/IC50 비교로, CRISPR arm은 각 오가노이드 자신의 비표적 대조 대비 log2FC로 평가하며, <strong>두 유전형 모두에서 죽는 유전자는 공통 필수</strong>로 분류되어 합성치사 근거가 되지 못합니다.`
      : `파이프라인의 예측 단계는 완료됐고, 이 단계는 실험 데이터를 기다립니다. 측정값이 들어오면 동일한 리포트에 자동으로 채워집니다.`,
    nodes);
}

/* Stage 6 - concordance */
{
  const nodes = [];
  if (concordance.length) {
    const hits = concordance.filter(c => c.experimental_hit === true || c.experimental_hit === "True").length;
    nodes.push(html(`<div class="callout"><div class="h">채점 결과</div>
      <p>공유 유전자 ${concordance.length}개 중 실험 적중 ${hits}개. 반증된 예측이 먼저 보고되며, 상위 집합이 전체 검정 집합과 같으면 enrichment는 검정 불가로 표시됩니다.</p></div>`));
    nodes.push(html(table([
      { title: "Gene", cls: "gene", render: c => esc(c.gene) },
      { title: "Predicted rank", render: c => c.predicted_rank },
      { title: "Predicted score", render: c => fmt(c.predicted_score, 3) },
      { title: "Measured", render: c => fmt(c.measured_effect, 3) },
      { title: "Hit", render: c => (c.experimental_hit === true || c.experimental_hit === "True") ? "✓" : "✗" }
    ], concordance)));
  } else {
    nodes.push(html(`<div class="callout"><div class="h">아직 채점 불가</div>
      <p>예측과 실험을 맞대는 단계입니다. <code>compare_prediction_with_experiment</code>가 반증된 예측 → 확인된 예측 → 예측이 놓친 적중 순으로 보고하고,
      precision@k를 permutation null에 대해 검정한 뒤 PREDICTIVE / WEAK / NOT PREDICTIVE 판정을 냅니다. 이 판정이 연구 질문에 대한 답입니다.</p></div>`));
  }
  stage("06", concordance.length ? "done" : "pending", "예측 vs 실험 — 연구 질문에 대한 답", 
    concordance.length ? "예측 순위와 측정값의 일치도입니다." : "실험 데이터가 들어오면 이 섹션이 연구 질문에 직접 답합니다.", nodes);
}

document.getElementById("foot").innerHTML = `
  <div>데이터 출처 — DepMap CRISPR gene effect (Chronos) &amp; model annotation, cBioPortal CCLE 변이/CN 콜,
  Broad Repurposing Hub, SynLethDB 2.0, NCBI PubMed (Entrez E-utilities).</div>
  <div style="margin-top:8px">한계 — 세포주 코호트는 환자 집단이 아니며 2D 배양에 적응한 모델에 편향돼 있습니다.
  단일 유전자 KO 의존성은 합성치사의 대리 지표이고, 세포주 수준의 선택성은 치료 창을 보장하지 않습니다.
  문헌 점수는 보고된 양이지 진실의 양이 아닙니다.</div>
  <div style="margin-top:8px">생성 ${esc(REPORT.generated_at)} · Biomni PDAC translational pipeline</div>`;
</script>
"""
