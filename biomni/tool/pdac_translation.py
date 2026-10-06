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
import urllib.parse
from datetime import datetime

from biomni.tool.synthetic_lethality import (
    _STRATIFY_ALIASES,
    _UTC,
    CBIOPORTAL_API,
    DEPLETION_THRESHOLD,
    PAN_ESSENTIAL_FRACTION,
    _annotate_copy_number_status,
    _annotate_mutation_status,
    _http_get,
    _http_post,
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

# The Repurposing Hub snapshot in the data lake predates the KRAS inhibitor era entirely: it has no
# sotorasib, no adagrasib, no MRTX compound and not a single SOS1 annotation. For a KRAS project that
# is not a gap in coverage, it is the whole pharmacology. ChEMBL is queried live to fill it.
# Curated RAS-pathway pharmacology, because neither local source carries it: the Repurposing Hub
# snapshot predates every KRAS inhibitor, and ChEMBL is a live API that is sometimes unavailable.
# Status is deliberately coarse ("approved" / "clinical" / "tool") rather than a precise trial phase,
# which moves faster than any snapshot can track. VERIFY BEFORE CITING - this is a starting point for
# a search, not a regulatory record.
CURATED_AS_OF = "2026-10"
CURATED_RAS_DRUGS = {
    "KRAS": [
        ("sotorasib (AMG 510)", "approved", "KRAS G12C covalent inhibitor"),
        ("adagrasib (MRTX849)", "approved", "KRAS G12C covalent inhibitor"),
        ("divarasib (GDC-6036)", "clinical", "KRAS G12C covalent inhibitor"),
        ("opnurasib (JDQ443)", "clinical", "KRAS G12C covalent inhibitor"),
        ("MRTX1133", "clinical", "KRAS G12D non-covalent inhibitor"),
        ("daraxonrasib (RMC-6236)", "clinical", "pan-RAS(ON) multi-selective inhibitor"),
        ("RMC-6291", "clinical", "KRAS G12C(ON) tri-complex inhibitor"),
        ("RMC-9805", "clinical", "KRAS G12D(ON) tri-complex inhibitor"),
    ],
    "SOS1": [
        ("BI-1701963", "clinical", "SOS1::KRAS interaction inhibitor"),
        ("MRTX0902", "clinical", "SOS1::KRAS interaction inhibitor"),
        ("BAY-293", "tool", "SOS1 inhibitor, chemical probe"),
    ],
    "PTPN11": [
        ("TNO155", "clinical", "SHP2 allosteric inhibitor"),
        ("RMC-4630", "clinical", "SHP2 allosteric inhibitor"),
        ("JAB-3068", "clinical", "SHP2 allosteric inhibitor"),
        ("SHP099", "tool", "SHP2 allosteric inhibitor, chemical probe"),
    ],
    "WRN": [
        ("HRO761", "clinical", "WRN helicase inhibitor (MSI-H context)"),
        ("VVD-133214", "clinical", "WRN helicase covalent inhibitor"),
    ],
    "PRMT5": [
        ("MRTX1719", "clinical", "MTA-cooperative PRMT5 inhibitor (MTAP-deleted context)"),
        ("AMG 193", "clinical", "MTA-cooperative PRMT5 inhibitor (MTAP-deleted context)"),
    ],
    "MAP2K1": [("trametinib", "approved", "MEK1/2 inhibitor"), ("selumetinib", "approved", "MEK1/2 inhibitor")],
    "MAPK1": [("ulixertinib", "clinical", "ERK1/2 inhibitor")],
    "CDK4": [("palbociclib", "approved", "CDK4/6 inhibitor"), ("abemaciclib", "approved", "CDK4/6 inhibitor"),
             ("ribociclib", "approved", "CDK4/6 inhibitor")],
    "CDK6": [("palbociclib", "approved", "CDK4/6 inhibitor"), ("abemaciclib", "approved", "CDK4/6 inhibitor")],
    "PTK2": [("defactinib", "clinical", "FAK inhibitor")],
    "EGFR": [("erlotinib", "approved", "EGFR inhibitor"), ("cetuximab", "approved", "EGFR monoclonal antibody")],
    "PARP1": [("olaparib", "approved", "PARP inhibitor"), ("niraparib", "approved", "PARP inhibitor")],
    "ATR": [("ceralasertib", "clinical", "ATR inhibitor"), ("berzosertib", "clinical", "ATR inhibitor")],
    "WEE1": [("adavosertib", "clinical", "WEE1 inhibitor")],
    "PKMYT1": [("lunresertib (RP-6306)", "clinical", "PKMYT1 inhibitor")],
    "YAP1": [("VT3989", "clinical", "TEAD palmitoylation inhibitor (YAP/TEAD axis)")],
    "TEAD1": [("VT3989", "clinical", "TEAD palmitoylation inhibitor"),
              ("IK-930", "clinical", "TEAD1-selective inhibitor"),
              ("K-975", "tool", "TEAD inhibitor, chemical probe")],
}
CURATED_STATUS_PHASE = {"approved": "Launched", "clinical": "Phase 1/Phase 2", "tool": "Preclinical"}

CHEMBL_API = "https://www.ebi.ac.uk/chembl/api/data"
CHEMBL_PHASE_NAMES = {4: "Launched", 3: "Phase 3", 2: "Phase 2", 1: "Phase 1", 0: "Preclinical"}

_DRUG_TABLE_CACHE: dict = {}
_CHEMBL_CACHE: dict = {}


def _chembl_targets_for_gene(gene: str) -> list | None:
    """ChEMBL target ids whose component really carries this gene symbol.

    The synonym search alone returns paralogs and complexes, so the gene symbol is re-checked
    against the component's GENE_SYMBOL synonyms before a target is accepted.
    """
    payload = _http_get(
        f"{CHEMBL_API}/target.json",
        params={
            "target_components__target_component_synonyms__component_synonym": gene,
            "limit": 20,
        },
        timeout=20,
        retries=1,  # a live API outage must not stall the whole candidate list
    )
    if payload is None:
        return None  # API failure, not an empty result - the caller must be able to tell them apart
    accepted = []
    for target in payload.get("targets", []):
        if target.get("organism") != "Homo sapiens":
            continue
        if target.get("target_type") != "SINGLE PROTEIN":
            continue  # complexes and PPI entries would attribute a drug to the wrong gene
        for component in target.get("target_components", []):
            symbols = {
                str(synonym.get("component_synonym", "")).upper()
                for synonym in component.get("target_component_synonyms", [])
                if synonym.get("syn_type") == "GENE_SYMBOL"
            }
            if gene.upper() in symbols:
                accepted.append(target["target_chembl_id"])
                break
    return accepted


def _drugs_from_chembl(gene: str) -> list | None:
    """Live ChEMBL lookup: compounds with a recorded mechanism of action against ``gene``.

    Returns None when the API cannot be reached, so an outage is never reported as "no drug exists".
    """
    gene = gene.strip().upper()
    if gene in _CHEMBL_CACHE:
        return _CHEMBL_CACHE[gene]

    targets = _chembl_targets_for_gene(gene)
    if targets is None:
        return None
    records = []
    for target_id in targets:
        mechanisms = _http_get(
            f"{CHEMBL_API}/mechanism.json",
            params={"target_chembl_id": target_id, "limit": 100},
            timeout=20,
            retries=1,
        )
        if mechanisms is None:
            # A target we know exists but whose mechanism list we could not read: reporting the
            # remaining genes as complete would turn an outage into "no compound exists".
            return None
        if not mechanisms.get("mechanisms"):
            continue
        by_molecule = {}
        for mechanism in mechanisms.get("mechanisms", []):
            molecule_id = mechanism.get("molecule_chembl_id")
            if not molecule_id:
                continue
            by_molecule[molecule_id] = {
                "moa": mechanism.get("mechanism_of_action") or "",
                "action_type": (mechanism.get("action_type") or "").title(),
                "max_phase": mechanism.get("max_phase"),
            }
        if not by_molecule:
            continue

        names = _http_get(
            f"{CHEMBL_API}/molecule.json",
            params={
                "molecule_chembl_id__in": ",".join(list(by_molecule)[:50]),
                "limit": 50,
                "only": "molecule_chembl_id,pref_name,max_phase",
            },
            timeout=20,
            retries=1,
        )
        if names is None:
            return None
        name_by_id = {}
        for molecule in names.get("molecules", []):
            name_by_id[molecule["molecule_chembl_id"]] = (
                molecule.get("pref_name") or molecule["molecule_chembl_id"],
                molecule.get("max_phase"),
            )

        for molecule_id, detail in by_molecule.items():
            name, molecule_phase = name_by_id.get(molecule_id, (molecule_id, None))
            raw_phase = detail["max_phase"] if detail["max_phase"] is not None else molecule_phase
            try:
                phase_number = int(float(raw_phase)) if raw_phase is not None else 0
            except (TypeError, ValueError):
                phase_number = 0
            phase = CHEMBL_PHASE_NAMES.get(phase_number, "Preclinical")
            records.append(
                {
                    "gene": gene,
                    # Keep the name exactly as ChEMBL has it; when there is no preferred name the
                    # ChEMBL id stands in, and lowercasing it would make it unsearchable.
                    "drug": str(name),
                    "clinical_phase": phase,
                    "phase_rank": CLINICAL_PHASE_RANK.get(phase, 1),
                    "moa": detail["moa"],
                    "indication": "",
                    "disease_area": "",
                    "source": f"ChEMBL ({molecule_id})",
                }
            )
    _CHEMBL_CACHE[gene] = records
    return records


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
                        "source": "Repurposing Hub",
                    }
                )
    frame = pd.DataFrame(records)
    bundle = {"table": frame, "path": path, "n_drugs": table.shape[0], "n_genes": frame["gene"].nunique()}
    _DRUG_TABLE_CACHE[resolved] = bundle
    return bundle


def _drugs_from_curated(gene: str) -> list:
    """Curated RAS-pathway compounds for ``gene`` (offline, version-controlled, explicitly dated)."""
    records = []
    for name, status, moa in CURATED_RAS_DRUGS.get(gene.strip().upper(), []):
        phase = CURATED_STATUS_PHASE[status]
        records.append(
            {
                "gene": gene.strip().upper(),
                "drug": name,
                "clinical_phase": phase,
                "phase_rank": CLINICAL_PHASE_RANK.get(phase, 1),
                "moa": moa,
                "indication": "",
                "disease_area": "",
                "source": f"curated RAS set ({CURATED_AS_OF})",
            }
        )
    return records


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
    include_curated: bool = True,
    include_chembl: bool = True,
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
    include_curated : bool, optional
        Include the curated RAS-pathway drug set shipped with this module (default: True). The local
        Repurposing Hub snapshot predates every KRAS inhibitor, so without this a KRAS project sees
        no pharmacology at all.
    include_chembl : bool, optional
        Query ChEMBL live for compounds with a recorded mechanism against each gene (default: True).
        Failures degrade silently to the offline sources and are reported in the provenance block.
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
    if bundle is None and not include_curated and not include_chembl:
        return (
            "FAILURE: broad_repurposing_hub_phase_moa_target_info.parquet was not found in the data lake "
            "and both other sources are disabled."
        )

    minimum_rank = CLINICAL_PHASE_RANK.get(min_clinical_phase, 1)
    table = bundle["table"] if bundle is not None else pd.DataFrame(
        columns=["gene", "drug", "clinical_phase", "phase_rank", "moa", "indication", "disease_area", "source"]
    )
    chembl_failed = []
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
        extra = []
        if include_curated:
            extra.extend(_drugs_from_curated(gene))
        if include_chembl:
            found = _drugs_from_chembl(gene)
            if found is None:
                chembl_failed.append(gene)
            else:
                extra.extend(found)
        combined = pd.concat([table, pd.DataFrame(extra)], ignore_index=True) if extra else table

        hits = combined[(combined["gene"] == gene) & (combined["phase_rank"] >= minimum_rank)].copy()
        hits["inhibitory"] = hits["moa"].map(_is_inhibitory)
        # Inhibitory mechanisms first, then clinical phase: a launched substrate is not a target drug.
        hits = hits.sort_values(["inhibitory", "phase_rank"], ascending=[False, False]).drop_duplicates("drug")
        synlethdb_note = ""
        if gene in partners:
            synlethdb_note = f"SynLethDB partner of {driver_gene.upper()} ({', '.join(partners[gene]['sources'])})"

        if len(hits) == 0:
            undruggable.append(gene)
            rows.append({"gene": gene, "drug": "", "clinical_phase": "", "moa": "", "inhibitory": False,
                         "source": "", "verdict": "CRISPR-ONLY", "synlethdb": synlethdb_note})
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
                         "source": hit.get("source", "Repurposing Hub"),
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
                origin = str(row.get("source", "")).split(" (")[0]
                log.append(
                    f"             - {row['drug']} [{row['clinical_phase']}] {row['moa']}{marker}  <{origin}>"
                )
        if verdict == "CRISPR-ONLY":
            searched = ["Repurposing Hub"]
            if include_curated:
                searched.append("curated RAS set")
            if include_chembl:
                searched.append("ChEMBL" + (" (UNAVAILABLE - may be incomplete)" if gene in chembl_failed else ""))
            log.append(f"             - no compound found in: {', '.join(searched)}")

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
    if chembl_failed:
        log.append(
            f"  - ChEMBL was unreachable for {len(chembl_failed)} gene(s) ({', '.join(chembl_failed[:8])}). "
            "A CRISPR-ONLY verdict for those genes means 'nothing in the offline sources', not 'no compound "
            "exists'. Re-run when the API is back before designing around it."
        )
    log.append(
        "  - The curated RAS set is a dated snapshot maintained in this repository, not a registry. "
        f"Status labels (approved / clinical / tool) are as of {CURATED_AS_OF} and must be re-checked against "
        "ClinicalTrials.gov or the drug label before they are cited."
    )
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
    if bundle is not None:
        log.append(
            f"  - Broad Repurposing Hub: {bundle['path']} ({bundle['n_drugs']} compounds, "
            f"{bundle['n_genes']} annotated target genes). This snapshot predates the KRAS inhibitor era - "
            "it contains no sotorasib, adagrasib or MRTX compound and no SOS1 annotation."
        )
    if include_curated:
        log.append(
            f"  - Curated RAS-pathway set bundled with this module, as of {CURATED_AS_OF}: "
            f"{sum(len(v) for v in CURATED_RAS_DRUGS.values())} compounds across "
            f"{len(CURATED_RAS_DRUGS)} genes. Status is coarse (approved / clinical / tool) and must be "
            "re-checked against ClinicalTrials.gov or the label before it is cited."
        )
    if include_chembl:
        if chembl_failed:
            log.append(
                f"  - ChEMBL live lookup: UNAVAILABLE for {', '.join(chembl_failed)} (API error). Those genes "
                "were scored from the offline sources only, so their drug list may be incomplete."
            )
        else:
            log.append(f"  - ChEMBL live lookup ({CHEMBL_API}): mechanism-of-action records per gene")
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

    from biomni.tool.synthetic_lethality import (
        _entrez_efetch,
        _entrez_esearch,
        _score_literature,
        sl_literature_queries,
    )

    evidence = {}
    delay = 0.4 if api_key else 0.75
    for gene in genes:
        queries = sl_literature_queries(gene, driver, disease)

        pmids: list[str] = []
        used = []
        for label, query in queries:
            hits = _entrez_esearch(query, max_papers, email, api_key)
            time.sleep(delay)
            new = [pmid for pmid in hits if pmid not in pmids]
            used.append(f"[{label}] {len(hits)} hit(s)")
            pmids.extend(new)  # every tier runs, as in validate_sl_candidates_with_pubmed

        # Same depth as validate_sl_candidates_with_pubmed, so the two report the same score.
        records = _entrez_efetch(pmids[: max_papers * 2], email, api_key)
        time.sleep(delay)
        scored = _score_literature(records, gene, driver, disease)

        # A PubMed hit for "KRAS AND <gene> AND pancreatic" is often a panel paper that merely lists
        # the gene. Count the records whose text really carries BOTH symbols, so a high volume score
        # built on co-mention cannot pass for evidence about the interaction.
        def mentions_pair(record, _gene=gene):
            text = f"{record.get('title', '')} {record.get('abstract', '')}".upper()
            return driver.upper() in text and _gene.upper() in text

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


# Why each gene could plausibly sit downstream of KRAS. CURATED PRIOR KNOWLEDGE as of the date
# below, used only to annotate - it is not evidence from this screen, and a gene without an entry
# is not thereby a worse candidate.
KRAS_MECHANISM_AS_OF = "2026-10"
KRAS_MECHANISM_NOTES = {
    "RAB10": "Rab GTPase of the secretory/recycling route. KRAS-mutant PDAC scavenges nutrients by "
    "macropinocytosis, a membrane-traffic-dependent process, which is the plausible route by which a "
    "trafficking GTPase becomes limiting in this genotype. No paper tests KRAS-RAB10 directly.",
    "DOCK5": "Guanine nucleotide exchange factor for RAC1. KRAS signals to RAC1 for actin remodelling "
    "and transformation, so a RAC-GEF is mechanistically on-pathway. Entirely untested in the literature.",
    "PIK3CA": "Catalytic subunit of PI3K, the parallel effector arm of RAS alongside RAF-MEK-ERK. When "
    "the MAPK arm is shut down by a KRAS inhibitor, PI3K-AKT sustains survival signalling - the documented "
    "resistance mechanism to oncogenic KRAS inhibition in PDAC.",
    "CTNNB1": "Beta-catenin, effector of WNT signalling. WNT and RAS pathways converge on MYC and cell-cycle "
    "entry in gastrointestinal cancers, and PDAC lines retain WNT dependence, but the link to KRAS here is "
    "transcriptional convergence rather than a direct signalling relay.",
    "KLF5": "Transcription factor driving proliferative and basal-type programmes in pancreatic epithelium, "
    "induced downstream of RAS-MAPK. Reported as required for KRAS-driven transformation, though at least one "
    "report contests the interaction.",
}


def _protocol_cascade(frame, q_threshold, min_cohens_d, min_pct_mutant, max_pct_all):
    """Apply the filter cascade step by step and return (steps, surviving frame).

    Each step reports how many genes remain, so the report can show where a protocol loses its
    candidates instead of only showing what survived.
    """
    steps, mask = [], frame["gene"].notna()
    for label, condition in (
        (f"q < {q_threshold}", frame["q_value"] < q_threshold),
        (f"Cohen's d < {min_cohens_d}", frame["cohens_d"] < min_cohens_d),
        (f"%mutant dependent > {min_pct_mutant}%", frame["pct_a_dependent"] > min_pct_mutant),
        (f"%all lines dependent < {max_pct_all}%", frame["pct_all_dependent"] < max_pct_all),
        ("not pan-essential", ~frame["pan_essential"]),
    ):
        mask = mask & condition
        steps.append({"label": label, "remaining": int(mask.sum())})
    return steps, frame[mask].sort_values("effect_difference")


def _experiment_panel(frame, candidates: list, n_negative: int = 5, n_assay: int = 3):
    """Build the knockout panel: candidates, predicted-negative controls, assay positive controls.

    The negative arm is drawn from genes the same scan predicts to be non-selective but that still
    have a measurable knockout effect - without them precision@k equals the base rate and the
    experiment cannot distinguish a real prediction from a lucky one.
    """
    negatives = frame[
        (frame["effect_difference"].abs() < 0.005)
        & (frame["p_one_sided"].between(0.45, 0.55))
        & (frame["pct_all_dependent"].between(10, 45))
        & (~frame["pan_essential"])
    ].nlargest(n_negative, "pct_all_dependent")
    assay = frame[(frame["pan_essential"]) & (frame["pct_all_dependent"] > 95)].nsmallest(n_assay, "mean_a")

    rows = [{"role": "candidate", "gene": gene} for gene in candidates]
    rows += [
        {"role": "negative control", "gene": row["gene"],
         "note": f"predicted non-selective (delta={row['effect_difference']:+.3f}, p={row['p_one_sided']:.2f}); "
                 f"{row['pct_all_dependent']:.0f}% of all lines depend on it, so the knockout is still measurable"}
        for _, row in negatives.iterrows()
    ]
    rows += [
        {"role": "assay control", "gene": row["gene"],
         "note": f"common essential ({row['pct_all_dependent']:.0f}% of all lines); confirms delivery and editing"}
        for _, row in assay.iterrows()
    ]
    rows.append({"role": "non-targeting", "gene": "NTC x3", "note": "normalisation baseline for every organoid"})
    return rows


def _druggable_near_misses(frame, passing_genes: set, q_threshold, min_cohens_d, data_lake_path, limit: int = 5):
    """Genes that miss the strict gates but have a clinical-stage inhibitor.

    A protocol that gates hard on statistics and then demands a drug can end up with candidates that
    have no chemical matter at all. Surfacing the genes just under the line that DO have one makes
    that trade-off explicit instead of leaving it to be discovered after the fact.
    """
    hub = _load_repurposing_hub(data_lake_path)
    table = hub["table"] if hub is not None else None
    relaxed = frame[
        (frame["q_value"] < 0.25)
        & (frame["cohens_d"] < -0.5)
        & (frame["pct_a_dependent"] > 30)
        & (frame["pct_all_dependent"] < 50)
        & (~frame["pan_essential"])
        & (~frame["gene"].isin(passing_genes))
    ].sort_values("effect_difference")

    rows = []
    for _, row in relaxed.iterrows():
        gene = row["gene"]
        drugs = _drugs_from_curated(gene)
        if table is not None:
            drugs += table[table["gene"] == gene].to_dict("records")
        inhibitors = [d for d in drugs if _is_inhibitory(d.get("moa", ""))]
        clinical = [d for d in inhibitors if d["phase_rank"] >= 2]
        if not clinical:
            continue
        pathway, rationale = _combination_pathway(gene)
        best = max(clinical, key=lambda d: d["phase_rank"])
        rows.append(
            {
                "gene": gene,
                "effect_difference": float(row["effect_difference"]),
                "cohens_d": float(row["cohens_d"]),
                "q_value": float(row["q_value"]),
                "pct_a_dependent": float(row["pct_a_dependent"]),
                "pct_all_dependent": float(row["pct_all_dependent"]),
                "missed_gate": ", ".join(
                    filter(None, [
                        f"q={row['q_value']:.3f} (needs <{q_threshold})" if row["q_value"] >= q_threshold else "",
                        f"d={row['cohens_d']:.2f} (needs <{min_cohens_d})" if row["cohens_d"] >= min_cohens_d else "",
                    ])
                ),
                "n_inhibitors": len(inhibitors),
                "n_clinical": len(clinical),
                "best_drug": best["drug"],
                "best_phase": best["clinical_phase"],
                "combination_pathway": pathway,
                "combination_rationale": rationale,
            }
        )
        if len(rows) >= limit:
            break
    return rows


def _combination_candidates(frame, literature: dict, data_lake_path, q_max=0.15, d_max=-0.5, lit_min=70):
    """Section 3: near-miss genes with real literature support, with their drugs and refuting PMIDs.

    These miss the strict statistical gates but are documented, which makes them combination
    hypotheses rather than discoveries: the evidence they rest on is other people's, not this screen's.
    """
    hub = _load_repurposing_hub(data_lake_path)
    table = hub["table"] if hub is not None else None
    rows = []
    for _, row in frame.sort_values("effect_difference").iterrows():
        gene = row["gene"]
        entry = literature.get(gene)
        if entry is None or entry["score"] < lit_min:
            continue
        drugs = _drugs_from_curated(gene)
        if table is not None:
            drugs += table[table["gene"] == gene].to_dict("records")
        inhibitors = [d for d in drugs if _is_inhibitory(d.get("moa", ""))]
        clinical = [d for d in inhibitors if d["phase_rank"] >= 2]
        pathway, _ = _combination_pathway(gene)
        rows.append({
            "gene": gene,
            "effect_difference": float(row["effect_difference"]),
            "cohens_d": float(row["cohens_d"]),
            "q_value": float(row["q_value"]),
            "literature_score": entry["score"],
            "refuting_pmids": entry.get("refuting_pmids", []),
            "n_clinical": len(clinical),
            "best_drug": max(clinical, key=lambda d: d["phase_rank"])["drug"] if clinical else "",
            "best_phase": max(clinical, key=lambda d: d["phase_rank"])["clinical_phase"] if clinical else "",
            "combination_pathway": pathway,
        })
    return rows


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
    protocol_csv: str | None = None,
    protocol_q_threshold: float = 0.05,
    protocol_min_cohens_d: float = -0.8,
    protocol_min_pct_mutant: float = 30.0,
    protocol_max_pct_all: float = 50.0,
    biovalidation_csv: str | None = None,
    combination_max_genes: int = 30,
    data_lake_path: str | None = None,
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
    protocol_csv : str, optional
        Scan table (e.g. pancancer_g12d.csv) to run the strict filter cascade over. Defaults to
        ``<run_dir>/pancancer_g12d.csv`` when it exists. Adds the protocol funnel, the final
        candidate list and the organoid knockout panel to the report.
    protocol_q_threshold, protocol_min_cohens_d, protocol_min_pct_mutant, protocol_max_pct_all : optional
        The cascade's gates (defaults: q < 0.05, Cohen's d < -0.8, >30% of the mutant arm dependent,
        <50% of all lines dependent).
    biovalidation_csv : str, optional
        Table from ``validate_candidates_biologically`` (default: <run_dir>/biovalidation.csv). Adds the
        PDAC-restricted effect, KRAS axis, TCGA survival, drug status and structure confidence per candidate.
    combination_max_genes : int, optional
        How many near-miss genes to query PubMed for when building the combination-candidate table
        (default: 30).
    data_lake_path : str, optional
        Directory holding the Repurposing Hub table, used to find druggable near-miss candidates.
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

    protocol_path = protocol_csv or os.path.join(run_dir, "pancancer_g12d.csv")
    protocol_frame = _read_optional_csv(protocol_path)
    if protocol_frame is not None and "cohens_d" in protocol_frame.columns:
        steps, passing = _protocol_cascade(
            protocol_frame, protocol_q_threshold, protocol_min_cohens_d,
            protocol_min_pct_mutant, protocol_max_pct_all,
        )
        # The driver is the positive control, not a candidate: finding that KRAS-mutant lines depend on
        # KRAS validates the pipeline and says nothing about a partner.
        all_passing = passing["gene"].tolist()
        control_gene = driver.upper() if driver.upper() in all_passing else None
        protocol_genes = [g for g in all_passing if g != control_gene]
        payload["protocol"] = {
            "source": os.path.basename(protocol_path),
            "n_tested": int(len(protocol_frame)),
            "gates": {
                "q": protocol_q_threshold, "cohens_d": protocol_min_cohens_d,
                "pct_mutant": protocol_min_pct_mutant, "pct_all": protocol_max_pct_all,
            },
            "cascade": steps,
            "control_gene": control_gene,
            "control_row": passing[passing["gene"] == control_gene].to_dict("records")[0] if control_gene else None,
            "passing": passing[passing["gene"] != control_gene].head(10).to_dict("records"),
            "panel": _experiment_panel(protocol_frame, protocol_genes[:5]),
            "mechanisms": {g: KRAS_MECHANISM_NOTES.get(g, "") for g in protocol_genes[:10]},
            "mechanism_as_of": KRAS_MECHANISM_AS_OF,
            "near_miss": _druggable_near_misses(
                protocol_frame, set(all_passing), protocol_q_threshold, protocol_min_cohens_d, data_lake_path
            ),
        }
        payload["sections"].append("protocol")
        # Literature is fetched for the protocol's own candidates and the druggable near-misses,
        # since those are the genes a reader has to decide between.
        near_genes = [row["gene"] for row in payload["protocol"]["near_miss"]]
        combination_pool = protocol_frame[
            (protocol_frame["q_value"] < 0.15)
            & (protocol_frame["cohens_d"] < -0.5)
            & (~protocol_frame["pan_essential"])
            & (protocol_frame["gene"] != driver.upper())
        ].nsmallest(combination_max_genes, "effect_difference")
        payload["protocol"]["combination_pool"] = combination_pool["gene"].tolist()
        candidates = list(dict.fromkeys(
            protocol_genes[:10] + near_genes + combination_pool["gene"].tolist() + candidates
        ))

    literature = {}
    if include_literature and candidates:
        literature = _literature_evidence(
            disease, driver, candidates, max_papers=max_papers_per_gene, email=email, api_key=api_key
        )
        payload["sections"].append("literature")
    payload["literature"] = literature

    if "protocol" in payload and literature:
        pool = protocol_frame[protocol_frame["gene"].isin(payload["protocol"].get("combination_pool", []))]
        payload["protocol"]["combination"] = _combination_candidates(pool, literature, data_lake_path)

    biovalidation = _read_optional_csv(biovalidation_csv or os.path.join(run_dir, "biovalidation.csv"))
    if biovalidation is not None:
        payload["biovalidation"] = biovalidation.to_dict("records")
        payload["sections"].append("biovalidation")

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
  --step-1:#86b6ef; --step-2:#6da7ec; --step-3:#5598e7; --step-4:#3987e5; --step-5:#2a78d6;
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
    --step-1:#1c5cab; --step-2:#256abf; --step-3:#2a78d6; --step-4:#3987e5; --step-5:#5598e7;
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
  --step-1:#1c5cab; --step-2:#256abf; --step-3:#2a78d6; --step-4:#3987e5; --step-5:#5598e7;
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
.cards{display:grid;gap:12px;margin-bottom:18px}
.card{background:var(--surface);border:1px solid var(--rule);border-radius:3px;padding:16px 18px;box-shadow:var(--shadow)}
.card header{display:flex;flex-wrap:wrap;align-items:baseline;gap:10px;margin-bottom:10px}
.card .g{font-family:var(--sans);font-weight:700;font-size:19px;letter-spacing:-.01em}
.card .tier{font-family:var(--sans);font-size:10.5px;font-weight:700;letter-spacing:.09em;text-transform:uppercase;
  padding:2px 7px;border-radius:2px;border:1px solid currentColor}
.card .tier.t1{color:var(--good)}
.evidence{display:grid;grid-template-columns:repeat(auto-fit,minmax(210px,1fr));gap:12px}
.ev{border-top:2px solid var(--rule);padding-top:8px}
.ev .h{font-family:var(--sans);font-size:10.5px;font-weight:700;letter-spacing:.09em;text-transform:uppercase;
  color:var(--muted);margin-bottom:4px}
.ev .v{font-family:var(--mono);font-size:12px;color:var(--ink-2);line-height:1.55}
.ev.ok{border-top-color:var(--good)} .ev.no{border-top-color:var(--axis)} .ev.warn{border-top-color:var(--crit)}
.mech{margin-top:12px;font-size:14.5px;color:var(--ink-2);max-width:68ch}
.mech b{color:var(--ink);font-weight:600}
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

/* ---------- chart: filter cascade funnel ---------- */
function funnel(steps, total) {
  const width = 760, rowHeight = 34, height = (steps.length + 1) * rowHeight + 8;
  const svg = svgEl("svg", { viewBox: `0 0 ${width} ${height}`, role: "img" });
  const labelWidth = 230, plot = width - labelWidth - 70;
  const rows = [{ label: "tested genes", remaining: total }].concat(steps);
  rows.forEach((row, index) => {
    const y = index * rowHeight;
    const w = Math.max(2, (row.remaining / total) * plot);
    const name = svgEl("text", { class: "lab", x: labelWidth - 12, y: y + 20, "text-anchor": "end" });
    name.textContent = row.label;
    svg.appendChild(name);
    svg.appendChild(svgEl("rect", {
      x: labelWidth, y: y + 7, width: w, height: 16, rx: 3,
      fill: index === 0 ? "var(--axis)" : `var(--step-${Math.min(index, 5)})`,
    }));
    const value = svgEl("text", { class: "val", x: labelWidth + w + 8, y: y + 20 });
    value.textContent = row.remaining.toLocaleString();
    svg.appendChild(value);
    const hit = svgEl("rect", { class: "hit", x: 0, y: y, width: width, height: rowHeight });
    const dropped = index === 0 ? 0 : rows[index - 1].remaining - row.remaining;
    hover(hit, `<b>${esc(row.label)}</b><br>남은 유전자 ${row.remaining.toLocaleString()}` +
      (index ? `<br>이 단계에서 탈락 ${dropped.toLocaleString()}` : ""));
    svg.appendChild(hit);
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

/* Protocol: strict filter cascade, final candidates, knockout panel */
const protocol = REPORT.protocol;
if (protocol) {
  const gates = protocol.gates;
  const nodes = [figure(funnel(protocol.cascade, protocol.n_tested),
    `${esc(protocol.source)} 의 ${protocol.n_tested.toLocaleString()}개 유전자에 필터를 순서대로 적용한 결과. ` +
    "막대에 마우스를 올리면 각 단계에서 몇 개가 탈락했는지 보입니다.", "")];

  const cards = document.createElement("div");
  cards.className = "cards";
  protocol.passing.slice(0, 5).forEach(row => {
    const gene = row.gene;
    const lit = (REPORT.literature || {})[gene];
    const drugRows = (REPORT.drugs || []).filter(d => d.gene === gene && d.drug &&
      (d.inhibitory === true || d.inhibitory === "True"));
    const clinical = drugRows.filter(d => !/preclinical/i.test(String(d.clinical_phase)));
    const refuting = lit && lit.refuting_pmids ? lit.refuting_pmids : [];
    const litClass = !lit ? "no" : (refuting.length ? "warn" : (lit.score >= 65 ? "ok" : "no"));
    const mech = (protocol.mechanisms || {})[gene] || "";
    const bio = (REPORT.biovalidation || []).find(b => b.gene === gene);
    const card = document.createElement("article");
    card.className = "card";
    card.innerHTML = `
      <header><span class="g">${esc(gene)}</span>
        <span class="tier t1">프로토콜 충족</span>
        ${refuting.length ? '<span class="tier" style="color:var(--crit)">반박 근거 있음</span>' : ""}</header>
      <div class="evidence">
        <div class="ev ok"><div class="h">통계 근거</div><div class="v">
          Δ ${fmt(row.effect_difference, 3)} · d ${fmt(row.cohens_d, 2)} · q ${sci(row.q_value)}<br>
          변이군 의존 ${fmt(row.pct_a_dependent, 0)}% · 전체 ${fmt(row.pct_all_dependent, 0)}%</div></div>
        <div class="ev ${litClass}"><div class="h">문헌 근거</div><div class="v">
          ${lit ? `${lit.score}/100 · 논문 ${lit.n_papers}편 (쌍 동시 ${lit.n_pair_papers ?? 0})<br>${esc(String(lit.interpretation).split(" - ")[0])}` : "조회 안 됨"}
          ${refuting.length ? `<br><b>반박 PMID:</b> ${refuting.map(p => `<a href="https://pubmed.ncbi.nlm.nih.gov/${esc(p)}/" target="_blank" rel="noopener">${esc(p)}</a>`).join(", ")}` : ""}
        </div></div>
        <div class="ev ${clinical.length ? "ok" : "no"}"><div class="h">약물 가용성</div><div class="v">
          ${drugRows.length ? `저해제 ${drugRows.length}개 (임상 ${clinical.length})<br>${esc(drugRows[0].drug)} [${esc(drugRows[0].clinical_phase)}]`
            : "저해 기전 화합물 없음 → CRISPR arm 전용"}</div></div>
        ${bio ? `
        <div class="ev ${Number(bio.pdac_d_vs_other) <= -0.5 ? "ok" : "no"}"><div class="h">PDAC 한정 효과</div><div class="v">
          야생형 대비 d ${fmt(bio.pdac_d_vs_wt, 2)}<br>기타 PDAC 대비 d ${fmt(bio.pdac_d_vs_other, 2)}
          ${Number(bio.pdac_d_vs_other) > -0.3 ? "<br><b>PDAC 내부에서 소실</b>" : ""}</div></div>
        <div class="ev ${bio.tcga_survival === "yes" ? "ok" : "no"}"><div class="h">TCGA PAAD 생존</div><div class="v">
          ${esc(String(bio.tcga_survival))} · log-rank p ${sci(bio.tcga_logrank_p)}</div></div>
        <div class="ev no"><div class="h">축 · 구조</div><div class="v">
          ${esc(String(bio.kras_axis || "-"))}<br>pLDDT ${fmt(bio.mean_plddt, 1)} (pocket 점수 아님)</div></div>` : ""}
      </div>
      ${mech ? `<p class="mech"><b>KRAS 연결 기전</b> — ${esc(mech)}</p>` : ""}`;
    cards.appendChild(card);
  });
  nodes.push(cards);

  if (protocol.control_row) {
    const c = protocol.control_row;
    nodes.push(html(`<div class="callout"><div class="h">양성 대조 — ${esc(protocol.control_gene)}</div>
      <p>드라이버 자신이 같은 필터를 통과했습니다 (Δ ${fmt(c.effect_difference, 3)}, d ${fmt(c.cohens_d, 2)},
      q ${sci(c.q_value)}, 변이군 ${fmt(c.pct_a_dependent, 0)}% 의존). 파이프라인이 작동한다는 증거이며,
      <b>후보 목록에서는 제외</b>합니다 — 드라이버는 파트너가 아닙니다.</p></div>`));
  }

  if ((protocol.near_miss || []).length) {
    const rows = protocol.near_miss.map(n => ({
      label: n.gene,
      value: Math.abs(n.effect_difference),
      valueLabel: `${n.n_clinical}개 임상약물`,
      color: "var(--series-2)",
      tip: `<b>${esc(n.gene)}</b><br>Δ ${fmt(n.effect_difference, 3)} · d ${fmt(n.cohens_d, 2)} · q ${sci(n.q_value)}<br>` +
           `미달 기준: ${esc(n.missed_gate)}<br>저해제 ${n.n_inhibitors}개 (임상 ${n.n_clinical})<br>` +
           `${esc(n.best_drug)} [${esc(n.best_phase)}]` + (n.combination_pathway ? `<br>병용 축: ${esc(n.combination_pathway)}` : ""),
    }));
    nodes.push(figure(barChart(rows, { tickDigits: 2, labelWidth: 100 }),
      "통계 기준은 못 넘었지만 <b>임상단계 저해제가 있는</b> 후보. 막대는 효과 차이의 절댓값, 우측은 임상 약물 수.", ""));
    nodes.push(html(`<div class="tablewrap"><table><thead><tr><th>Gene</th><th>Δ</th><th>d</th><th>q</th>
      <th>미달 기준</th><th>약물</th><th>병용 축</th></tr></thead><tbody>` +
      protocol.near_miss.map(n => `<tr><td class="gene">${esc(n.gene)}</td><td>${fmt(n.effect_difference,3)}</td>
        <td>${fmt(n.cohens_d,2)}</td><td>${sci(n.q_value)}</td><td>${esc(n.missed_gate)}</td>
        <td>${esc(n.best_drug)} [${esc(n.best_phase)}]</td><td>${esc(n.combination_pathway || "-")}</td></tr>`).join("") +
      `</tbody></table></div>`));
    nodes.push(html(`<div class="callout"><div class="h">프로토콜이 드러낸 긴장</div>
      <p>통계가 가장 강한 후보에는 화합물이 없고, 약물이 풍부한 후보는 통계 문턱 아래에 있습니다.
      선택지는 셋입니다 — ① 기준 유지(CRISPR arm만), ② 기준을 <b>데이터를 보기 전에</b> 사전 등록해 완화,
      ③ CRISPR arm과 약물 arm을 분리해 서로 다른 질문으로 검증. 셋째를 권합니다: 두 arm의 불일치 자체가 정보입니다.</p></div>`));
  }

  const body = stage("05", "done", "프로토콜 적용 — 최종 후보",
    `사전 지정한 필터(q &lt; ${gates.q}, Cohen's d &lt; ${gates.cohens_d}, 변이군 의존 &gt; ${gates.pct_mutant}%, ` +
    `전체 의존 &lt; ${gates.pct_all}%, pan-essential 제외)를 ${protocol.n_tested.toLocaleString()}개 유전자에 적용했습니다. ` +
    `<strong>통과 ${protocol.cascade[protocol.cascade.length - 1].remaining}개</strong>` +
    (protocol.control_gene ? ` (드라이버 ${esc(protocol.control_gene)} 양성대조 포함, 후보는 ${protocol.passing.length}개).` : ".") + " " +
    "각 후보마다 통계·문헌·약물·기전 네 가지 근거를 따로 제시하며, 어느 하나가 비어 있으면 비어 있는 대로 표시합니다.",
    nodes);
  body.appendChild(html(table([
    { title: "Gene", cls: "gene", render: r => esc(r.gene) },
    { title: "MUT mean", render: r => fmt(r.mean_a, 3) },
    { title: "WT mean", render: r => fmt(r.mean_b, 3) },
    { title: "Δ", render: r => fmt(r.effect_difference, 3) },
    { title: "Cohen's d", render: r => fmt(r.cohens_d, 2) },
    { title: "q", render: r => sci(r.q_value) },
    { title: "%MUT dep", render: r => fmt(r.pct_a_dependent, 0) },
    { title: "%all dep", render: r => fmt(r.pct_all_dependent, 0) },
  ], protocol.passing)));

  const combination = protocol.combination || [];
  if (combination.length) {
    const nodes3 = [html(`<div class="tablewrap"><table><thead><tr><th>Gene</th><th>Δ</th><th>d</th><th>q</th>
      <th>문헌</th><th>반박 PMID</th><th>임상약물</th><th>병용 축</th></tr></thead><tbody>` +
      combination.map(c => `<tr><td class="gene">${esc(c.gene)}</td><td>${fmt(c.effect_difference,3)}</td>
        <td>${fmt(c.cohens_d,2)}</td><td>${fmt(c.q_value,3)}</td><td>${c.literature_score}</td>
        <td>${c.refuting_pmids.length ? c.refuting_pmids.map(p => `<a href="https://pubmed.ncbi.nlm.nih.gov/${esc(p)}/" target="_blank" rel="noopener">${esc(p)}</a>`).join(" ") : "—"}</td>
        <td>${c.n_clinical ? `${c.n_clinical}개 · ${esc(c.best_drug)} [${esc(c.best_phase)}]` : "없음"}</td>
        <td>${esc(c.combination_pathway || "-")}</td></tr>`).join("") + `</tbody></table></div>`)];
    stage("06", "done", "병용 후보 — 통계 기준 미달이나 문헌이 뒷받침하는 유전자",
      "q &lt; 0.15, Cohen's d &lt; −0.5, 문헌 점수 ≥ 70을 모두 만족하는 유전자입니다. " +
      "이들이 기대는 근거는 <strong>이 스크린이 아니라 다른 연구자들의 것</strong>이므로 발견이 아니라 병용 가설로 읽어야 합니다. " +
      "반박 PMID가 있는 후보는 작동 전에 그 논문부터 확인하세요.", nodes3);
  }

  const noDrug = (REPORT.biovalidation || []).filter(b => !Number(b.n_clinical_inhibitors));
  if (noDrug.length) {
    const nodes4 = [html(`<div class="tablewrap"><table><thead><tr><th>후보</th><th>1-hop 약물 가능 이웃 (STRING ≥ 0.7)</th></tr></thead>
      <tbody>${noDrug.map(b => `<tr><td class="gene">${esc(b.gene)}</td>
        <td>${(b.druggable_neighbours && String(b.druggable_neighbours) !== "nan")
          ? esc(String(b.druggable_neighbours)) : "없음 — CRISPR·degrader 외 경로 없음"}</td></tr>`).join("")}
      </tbody></table></div>`)];
    nodes4.push(html(`<div class="callout"><div class="h">읽는 법</div>
      <p>이웃을 저해하는 것은 <b>다른 노드를 때리는 것</b>이라 후보 자체의 검증이 되지 않습니다.
      경로 가설로는 쓸 수 있지만, 그 후보가 필요한지에 대한 직접 답은 여전히 녹아웃입니다.</p></div>`));
    stage("07", "done", "약물 없는 후보 — 네트워크 이웃", "", nodes4);
  }

  const panel = protocol.panel || [];
  const roleLabel = { "candidate": "검증 대상", "negative control": "음성 대조", "assay control": "어세이 대조", "non-targeting": "비표적 대조" };
  const panelNodes = [html(`<div class="callout"><div class="h">사전 지정 기준 (측정 전 고정)</div>
    <p>적중: KRAS G12D 오가노이드에서 <b>log2FC ≤ −0.5</b> <em>그리고</em> 야생형 대비 선택성 <b>Δ ≤ −0.5</b>.
    정상 췌장 오가노이드에서 log2FC ≤ −0.5면 독성으로 탈락(치료 창 없음).
    음성 대조 5개 중 2개 이상이 적중 기준을 넘으면 어세이 무효.</p></div>`)];
  panelNodes.push(html(`<div class="tablewrap"><table><thead><tr><th>역할</th><th>유전자</th><th>근거</th></tr></thead>
    <tbody>${panel.map(p => `<tr><td>${esc(roleLabel[p.role] || p.role)}</td>
      <td class="gene">${esc(p.gene)}</td><td>${esc(p.note || "")}</td></tr>`).join("")}</tbody></table></div>`));
  panelNodes.push(html(`<div class="callout"><div class="h">Isogenic KRAS G12D — KI vs KO 3-arm</div>
    <p>세포주 패널은 유전형과 배경이 교란되므로, 동일 배경에서 KRAS만 바꾸는 축이 필요합니다.
    ① <b>KI arm</b>: 야생형 정상 췌장 오가노이드에 G12D knock-in → 파트너 KO 효과가 <em>생겨나는지</em>.
    ② <b>KO arm</b>: G12D PDO에서 allele 특이 sgRNA로 변이 대립유전자만 제거 → 효과가 <em>사라지는지</em>.
    ③ <b>약리학적 phenocopy arm</b>: 같은 PDO에 MRTX1133 처리 → 유전적 KO arm과 같은 방향이면 mutant KRAS 의존성이 이중 확인됩니다.
    세 arm이 일치해야 'G12D 의존적'이라 말할 수 있습니다 — KI만으로는 과발현 인공산물을, KO만으로는 적응 효과를 배제하지 못합니다.</p></div>`));
  stage("08", "done", "오가노이드 녹아웃 패널",
    "음성 대조를 <strong>예측이 낮게 평가한 유전자</strong>로 넣는 것이 핵심입니다. 상위 후보만 검증하면 " +
    "precision@k가 base rate와 같아져 '맞췄다'를 '우연보다 낫다'와 구분할 수 없습니다. " +
    "어세이 대조는 공통 필수 유전자로 전달·편집 효율을 확인합니다.", panelNodes);
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
  stage("09", crispr.length || curves.length ? "done" : "pending", "오가노이드 기능 검증",
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
  stage("10", concordance.length ? "done" : "pending", "예측 vs 실험 — 연구 질문에 대한 답",
    concordance.length ? "예측 순위와 측정값의 일치도입니다." : "실험 데이터가 들어오면 이 섹션이 연구 질문에 직접 답합니다.", nodes);
}

const gaps = [];
if (!(REPORT.biovalidation || []).length) gaps.push("생물학적 검증 미실행");
gaps.push("Sanger Project SCORE 독립 복제 미실행 — 데이터 레이크에 없어 현재 후보는 단일 스크린 유래입니다");
gaps.push("AlphaFold pocket quality 미산출 — pocket detector(fpocket 등)가 없어 pLDDT(모델 신뢰도)로 대체했습니다");
if ((REPORT.biovalidation || []).some(b => !b.chembl_class)) gaps.push("ChEMBL target class 일부 미확보 — API 장애 시점의 공백입니다");
gaps.push("야생형 대조군은 유전형만 야생형이며 NF1 결손·BRAF·RTK 증폭으로 RAS 경로가 켜진 세포주를 포함합니다");
gaps.push("TCGA 생존은 bulk 종양의 연관이라 기질 기여를 분리하지 못하며, 필요성의 증거가 아닙니다");
document.getElementById("stages").appendChild(html(`<section class="stage"><div class="rail">
  <span class="num">!</span><span class="state pending">한계</span></div>
  <div class="body"><h2>남은 공백</h2>
  <div class="finding">이 분석이 <strong>하지 못한</strong> 것들입니다. 결과를 읽을 때 함께 읽어야 합니다.</div>
  <ul class="refs">${gaps.map(g => `<li>${g}</li>`).join("")}</ul></div></section>`).firstChild);

document.getElementById("foot").innerHTML = `
  <div>데이터 출처 — DepMap CRISPR gene effect (Chronos) &amp; model annotation, cBioPortal CCLE 변이/CN 콜,
  Broad Repurposing Hub, SynLethDB 2.0, NCBI PubMed (Entrez E-utilities).</div>
  <div style="margin-top:8px">한계 — 세포주 코호트는 환자 집단이 아니며 2D 배양에 적응한 모델에 편향돼 있습니다.
  단일 유전자 KO 의존성은 합성치사의 대리 지표이고, 세포주 수준의 선택성은 치료 창을 보장하지 않습니다.
  문헌 점수는 보고된 양이지 진실의 양이 아닙니다.</div>
  <div style="margin-top:8px">생성 ${esc(REPORT.generated_at)} · Biomni PDAC translational pipeline</div>`;
</script>
"""


# ---------------------------------------------------------------------------
# Alternative stratifications: the KRAS mutant/wild-type dichotomy has no control arm
# ---------------------------------------------------------------------------
def _welch_scan(gene_effect, group_a_ids: list, group_b_ids: list, min_valid: int = 3):
    """Genome-wide one-sided Welch scan of group A vs group B on the gene-effect matrix.

    Returns a DataFrame with one row per gene: group means, effect difference (A - B, negative =
    A is more dependent), Cohen's d, one-sided p for "A more depleted", BH q, and dependent
    fractions. Shared by every stratification tool so the statistics cannot drift between them.
    """
    import numpy as np
    import pandas as pd
    from scipy import stats

    from biomni.tool.synthetic_lethality import _benjamini_hochberg

    block_a = gene_effect.loc[group_a_ids]
    block_b = gene_effect.loc[group_b_ids]
    usable = block_a.columns[(block_a.notna().sum() >= min_valid) & (block_b.notna().sum() >= min_valid)]
    block_a, block_b = block_a[usable], block_b[usable]

    tstat, p_two = stats.ttest_ind(block_a.values, block_b.values, axis=0, equal_var=False, nan_policy="omit")
    tstat = np.asarray(tstat, dtype=float)
    p_two = np.where(np.isfinite(np.asarray(p_two, dtype=float)), np.asarray(p_two, dtype=float), 1.0)
    p_one = np.where(tstat < 0, p_two / 2.0, 1.0 - p_two / 2.0)

    mean_a, mean_b = np.asarray(block_a.mean()), np.asarray(block_b.mean())
    pooled = np.sqrt((np.asarray(block_a.std(ddof=1)) ** 2 + np.asarray(block_b.std(ddof=1)) ** 2) / 2)
    with np.errstate(divide="ignore", invalid="ignore"):
        cohens_d = np.where(pooled > 0, (mean_a - mean_b) / pooled, np.nan)

    pan_fraction = (gene_effect[usable] < DEPLETION_THRESHOLD).sum() / gene_effect[usable].notna().sum()
    return pd.DataFrame(
        {
            "gene": [c.split(" (")[0] for c in usable],
            "n_a": np.asarray(block_a.notna().sum()),
            "n_b": np.asarray(block_b.notna().sum()),
            "mean_a": mean_a,
            "mean_b": mean_b,
            "effect_difference": mean_a - mean_b,
            "cohens_d": cohens_d,
            "p_one_sided": p_one,
            "q_value": _benjamini_hochberg(p_one),
            "pct_a_dependent": np.asarray((block_a < DEPLETION_THRESHOLD).sum() / block_a.notna().sum()) * 100,
            "pct_b_dependent": np.asarray((block_b < DEPLETION_THRESHOLD).sum() / block_b.notna().sum()) * 100,
            "pct_all_dependent": np.asarray(pan_fraction) * 100,
        }
    ).assign(pan_essential=lambda f: f["pct_all_dependent"] >= PAN_ESSENTIAL_FRACTION * 100)


def _apply_effect_filters(frame, min_cohens_d, min_pct_mutant_dependent, max_pct_all_dependent):
    """Optional effect-size and selectivity gates on a _welch_scan table.

    These are separate from the significance gates on purpose: a p-value says the difference is
    unlikely to be noise, these say it is large, it applies to most of the altered group, and it is
    not simply a gene that most cell lines need.
    """
    mask = frame["gene"].notna()
    if min_cohens_d is not None:
        mask &= frame["cohens_d"] <= min_cohens_d
    if min_pct_mutant_dependent is not None:
        mask &= frame["pct_a_dependent"] >= min_pct_mutant_dependent
    if max_pct_all_dependent is not None:
        mask &= frame["pct_all_dependent"] <= max_pct_all_dependent
    return mask


def _allele_groups(annotated, allele: str | None):
    """Split an annotated cohort into (allele-mutant, wild-type, other-allele) ModelID lists.

    Lines carrying a different allele of the same driver are excluded from both arms rather than
    pooled into either: they are neither the genotype under test nor a clean control.
    """
    mutant = annotated[annotated["MutationStatus"] == "MUT"]
    wildtype = annotated[annotated["MutationStatus"] == "WT"]
    if not allele:
        return mutant["ModelID"].tolist(), wildtype["ModelID"].tolist(), []
    carries = mutant["Variant"].astype(str).str.contains(allele, case=False, na=False)
    return (
        mutant[carries]["ModelID"].tolist(),
        wildtype["ModelID"].tolist(),
        mutant[~carries]["ModelID"].tolist(),
    )


def discover_sl_pan_cancer_with_context(
    driver_gene: str = "KRAS",
    allele: str | None = "G12D",
    context_cancer_type: str = "Pancreatic Cancer",
    data_lake_path: str | None = None,
    mutation_csv_path: str | None = None,
    top_n: int = 20,
    p_threshold: float = 0.05,
    fdr_threshold: float = 0.25,
    min_effect_difference: float = -0.2,
    max_mutant_mean_effect: float = -0.3,
    min_cohens_d: float | None = None,
    min_pct_mutant_dependent: float | None = None,
    max_pct_all_dependent: float | None = None,
    exclude_pan_essential: bool = True,
    output_csv_path: str | None = None,
) -> str:
    """Find driver-selective dependencies pan-cancer, then ask whether they hold in one cancer type.

    In PDAC the KRAS contrast has no control arm - 40 mutant lines against 4 wild-type - so a
    within-PDAC scan cannot clear a genome-wide FDR whatever the biology is. Pooling every lineage
    fixes the power problem (KRAS G12D: 47 mutant vs 791 wild-type) and creates a new one: a hit can
    be a lineage effect rather than a genotype effect. This tool therefore runs two steps that
    answer different questions, and reports both:

    1. **Pan-cancer selectivity** - is the dependency stronger in driver-mutant lines than in
       wild-type lines, across all lineages? Each candidate is re-tested with its dominant lineage
       removed, so a hit carried by one tissue is visible as such.
    2. **Context presence** - is the dependency actually there in ``context_cancer_type`` lines of
       that genotype? This is a magnitude question (how deep, in how many lines), not a second
       genotype contrast, because the context cohort cannot support one.

    Lines carrying a different allele of the same driver are excluded from both arms.

    Parameters
    ----------
    driver_gene : str, optional
        Driver defining the genotype, e.g. "KRAS" (default).
    allele : str, optional
        Restrict the mutant arm to one protein change, e.g. "G12D" (default). None uses any mutation.
    context_cancer_type : str, optional
        Cancer type the candidates are re-checked in (default: "Pancreatic Cancer").
    data_lake_path, mutation_csv_path : str, optional
        DepMap directory and an optional genotype table overriding the default call source.
    top_n : int, optional
        Number of candidates carried into the context step and printed (default: 20).
    p_threshold, fdr_threshold, min_effect_difference, max_mutant_mean_effect : optional
        Significance and effect-size gates for the pan-cancer step.
    min_cohens_d : float, optional
        Require Cohen's d at or below this value (e.g. -0.8 for a large effect). None disables it.
    min_pct_mutant_dependent : float, optional
        Require this percentage of the altered arm to be dependent (gene effect < the depletion
        threshold). Separates a shift of the whole group from a shift driven by a few lines.
    max_pct_all_dependent : float, optional
        Reject genes depleted in more than this percentage of ALL screened lines, a stricter and
        more explicit control than the pan-essential flag alone.
    exclude_pan_essential : bool, optional
        Drop common-essential genes (default: True).
    output_csv_path : str, optional
        Write the full pan-cancer table to this CSV path.

    Returns
    -------
    str
        A research log with both cohorts, the ranked candidates with their pan-cancer statistics,
        the leave-one-lineage-out check, the context verdict per gene, and a CANDIDATE_GENES line.

    """
    import numpy as np
    import pandas as pd

    driver = driver_gene.strip().upper()
    label = f"{driver} {allele}" if allele else f"{driver} mutant"
    log = [
        "=" * 78,
        f"PAN-CANCER DISCOVERY + {context_cancer_type.upper()} CONTEXT CHECK - {label}",
        f"Run at {datetime.now(tz=_UTC).strftime('%Y-%m-%d %H:%M UTC')}",
        "=" * 78,
        "",
    ]

    try:
        bundle = _load_depmap(data_lake_path)
    except FileNotFoundError as e:
        return f"FAILURE: {e}"
    gene_effect = bundle["gene_effect"]
    screened = bundle["model"][bundle["model"]["ModelID"].isin(gene_effect.index)]

    try:
        annotated, source = _annotate_mutation_status(screened, driver, bundle["data_lake_path"], mutation_csv_path)
    except (RuntimeError, ValueError) as e:
        return f"FAILURE: {e}"

    mutant_ids, wildtype_ids, other_allele_ids = _allele_groups(annotated, allele)
    if len(mutant_ids) < MIN_GROUP_SIZE or len(wildtype_ids) < MIN_GROUP_SIZE:
        return (
            f"FAILURE: pan-cancer groups too small ({len(mutant_ids)} vs {len(wildtype_ids)}). "
            f"Check the allele spelling ('{allele}') against the variant strings in the mutation source."
        )

    lineages = annotated[annotated["ModelID"].isin(mutant_ids)]["OncotreeLineage"].value_counts()
    log.append("STEP 1 | Pan-cancer cohort")
    log.append(f"  Genotype source: {source}")
    log.append(f"  {label}: {len(mutant_ids)} lines | {driver} wild-type: {len(wildtype_ids)} lines")
    if other_allele_ids:
        log.append(f"  Excluded from both arms (other {driver} allele): {len(other_allele_ids)} lines")
    log.append(f"  Mutant-arm lineages: {lineages.head(6).to_dict()}")

    results = _welch_scan(gene_effect, mutant_ids, wildtype_ids)
    selected = results[
        (results["p_one_sided"] < p_threshold)
        & (results["q_value"] < fdr_threshold)
        & (results["effect_difference"] <= min_effect_difference)
        & (results["mean_a"] <= max_mutant_mean_effect)
        & _apply_effect_filters(results, min_cohens_d, min_pct_mutant_dependent, max_pct_all_dependent)
    ]
    if exclude_pan_essential:
        selected = selected[~selected["pan_essential"]]
    selected = selected.sort_values(["effect_difference", "p_one_sided"]).head(top_n).reset_index(drop=True)

    log.append("")
    log.append("STEP 2 | Pan-cancer selectivity scan (one-sided Welch, BH-FDR)")
    log.append(f"  Genes tested: {len(results)}")
    log.append(f"  p < {p_threshold}: {(results['p_one_sided'] < p_threshold).sum()} | "
               f"q < {fdr_threshold}: {(results['q_value'] < fdr_threshold).sum()}")
    log.append(f"  Surviving all gates: {len(selected)}")

    # Leave-one-lineage-out: a candidate carried by a single tissue is a lineage effect wearing a
    # genotype label. Re-run the contrast with the mutant arm's dominant lineage removed.
    dominant = lineages.index[0] if len(lineages) else None
    lineage_delta = {}
    if dominant is not None and len(selected):
        keep_mutant = annotated[
            annotated["ModelID"].isin(mutant_ids) & (annotated["OncotreeLineage"] != dominant)
        ]["ModelID"].tolist()
        keep_wildtype = annotated[
            annotated["ModelID"].isin(wildtype_ids) & (annotated["OncotreeLineage"] != dominant)
        ]["ModelID"].tolist()
        if len(keep_mutant) >= MIN_GROUP_SIZE and len(keep_wildtype) >= MIN_GROUP_SIZE:
            columns = [bundle["gene_columns"][g] for g in selected["gene"] if g in bundle["gene_columns"]]
            subset = gene_effect[columns]
            delta = subset.loc[keep_mutant].mean() - subset.loc[keep_wildtype].mean()
            lineage_delta = {c.split(" (")[0]: float(v) for c, v in delta.items()}
            log.append(
                f"  Leave-one-lineage-out: dominant lineage '{dominant}' removed "
                f"({len(keep_mutant)} vs {len(keep_wildtype)} lines remain)"
            )

    # Context step: is the dependency present in this cancer type's lines of the same genotype?
    context_cohort, context_note = _select_cancer_models(bundle["model"], context_cancer_type)
    context_cohort = context_cohort[context_cohort["ModelID"].isin(gene_effect.index)]
    context_annotated, _ = _annotate_mutation_status(
        context_cohort, driver, bundle["data_lake_path"], mutation_csv_path
    )
    context_mutant, context_wildtype, _ = _allele_groups(context_annotated, allele)

    log.append("")
    log.append(f"STEP 3 | Context check in {context_cancer_type}")
    log.append(f"  {context_note}")
    log.append(f"  {label} lines in context: {len(context_mutant)} | {driver} wild-type: {len(context_wildtype)}")
    if len(context_mutant) < MIN_GROUP_SIZE:
        log.append(
            f"  Too few {label} lines in {context_cancer_type} to check presence; the pan-cancer result "
            "stands on its own and needs an orthogonal model."
        )

    rows = []
    for _, row in selected.iterrows():
        column = bundle["gene_columns"].get(row["gene"])
        context_mean, context_pct = np.nan, np.nan
        verdict = "NOT CHECKED"
        if column is not None and len(context_mutant) >= MIN_GROUP_SIZE:
            values = gene_effect.loc[context_mutant, column].dropna()
            if len(values):
                context_mean = float(values.mean())
                context_pct = float((values < DEPLETION_THRESHOLD).mean() * 100)
                if context_mean <= DEPLETION_THRESHOLD and context_pct >= 50:
                    verdict = "PRESENT IN CONTEXT"
                elif context_mean <= max_mutant_mean_effect:
                    verdict = "WEAK IN CONTEXT"
                else:
                    verdict = "ABSENT IN CONTEXT"
        rows.append(
            {
                **row.to_dict(),
                "delta_excl_dominant_lineage": lineage_delta.get(row["gene"], np.nan),
                "context_mean_effect": context_mean,
                "context_pct_dependent": context_pct,
                "context_verdict": verdict,
            }
        )
    table = pd.DataFrame(rows)

    log.append("")
    log.append(f"TOP {len(table)} CANDIDATES")
    log.append(
        f"{'gene':<12}{'MUT':>8}{'WT':>8}{'delta':>8}{'q':>9}{'no-lin':>8}"
        f"{'ctx mean':>10}{'ctx dep':>9}  context verdict"
    )
    log.append("-" * 92)
    for _, row in table.iterrows():
        no_lineage = "n/a" if not np.isfinite(row["delta_excl_dominant_lineage"]) else f"{row['delta_excl_dominant_lineage']:.3f}"
        context_mean = "n/a" if not np.isfinite(row["context_mean_effect"]) else f"{row['context_mean_effect']:.3f}"
        context_pct = "n/a" if not np.isfinite(row["context_pct_dependent"]) else f"{row['context_pct_dependent']:.0f}%"
        log.append(
            f"{row['gene']:<12}{row['mean_a']:>8.3f}{row['mean_b']:>8.3f}{row['effect_difference']:>8.3f}"
            f"{row['q_value']:>9.4f}{no_lineage:>8}{context_mean:>10}{context_pct:>9}  {row['context_verdict']}"
        )

    confirmed = table[table["context_verdict"] == "PRESENT IN CONTEXT"]["gene"].tolist()
    log.append("")
    log.append(f"CANDIDATE_GENES: {', '.join(confirmed) or ', '.join(table['gene'].head(top_n))}")
    log.append(
        f"  ({len(confirmed)} of {len(table)} candidates are actually dependent in {context_cancer_type} "
        f"{label} lines; the rest are pan-cancer genotype effects that this context does not show)"
    )

    log.append("")
    log.append("QC WARNINGS")
    log.append(
        f"  - Pooling lineages is what makes this test possible and is also its main threat: the mutant arm is "
        f"{lineages.iloc[0] / len(mutant_ids):.0%} {dominant}. The 'no-lin' column is the same contrast with that "
        "lineage removed - if it collapses toward zero, the hit was a tissue effect."
    )
    log.append(
        "  - The context column is a magnitude check, not a genotype contrast. It says the dependency exists in "
        f"those lines, not that it is selective for {driver} within {context_cancer_type}."
    )
    log.append(
        f"  - {len(other_allele_ids)} lines with a different {driver} allele were excluded from both arms. "
        "Pooling them into either arm would blur an allele-specific effect."
    )
    log.append(
        "  - Wild-type here means 'profiled and not mutated', which still includes lines with RAS-pathway "
        "activation by other means (NF1 loss, BRAF, RTK amplification)."
    )

    if output_csv_path:
        os.makedirs(os.path.dirname(os.path.abspath(output_csv_path)), exist_ok=True)
        results.sort_values(["effect_difference", "p_one_sided"]).to_csv(output_csv_path, index=False)
        log.append(f"  Full pan-cancer table written to {output_csv_path}")

    log.append("")
    log.append("PROVENANCE")
    for item in bundle["provenance"]:
        log.append(f"  - {item}")
    log.append(f"  - Genotype calls: {source}")
    log.append(
        f"  - Statistics: one-sided Welch t-test on Chronos gene effect, BH-FDR across {len(results)} genes; "
        f"depletion threshold {DEPLETION_THRESHOLD}, pan-essential cutoff {PAN_ESSENTIAL_FRACTION:.0%}"
    )
    return "\n".join(log)


def discover_comutation_stratified_sl(
    driver_gene: str = "KRAS",
    comutation_gene: str = "TP53",
    comutation_event: str = "mutation",
    driver_allele: str | None = None,
    cancer_type: str = "Pancreatic Cancer",
    data_lake_path: str | None = None,
    mutation_csv_path: str | None = None,
    top_n: int = 20,
    p_threshold: float = 0.05,
    fdr_threshold: float = 0.25,
    min_effect_difference: float = -0.2,
    max_mutant_mean_effect: float = -0.3,
    min_cohens_d: float | None = None,
    min_pct_mutant_dependent: float | None = None,
    max_pct_all_dependent: float | None = None,
    exclude_pan_essential: bool = True,
    output_csv_path: str | None = None,
) -> str:
    """Stratify driver-mutant lines by a second alteration, instead of by the driver itself.

    When a driver is near-universal in a cancer type its wild-type arm vanishes, and no amount of
    statistics recovers a contrast that the panel does not contain. The co-mutation design asks a
    different question that the same panel *can* answer: among lines that all carry the driver,
    does a second alteration change what they depend on? In PDAC that gives usable arms - KRAS-mutant
    lines split by TP53, SMAD4 or CDKN2A status - and the answer is directly actionable, because
    those co-alterations are what distinguishes one patient's tumour from another's.

    The result is a dependency of the *co-altered genotype*, not of the driver. The report states
    that explicitly: a hit here means "KRAS-mutant lines that also lost SMAD4 need this gene", which
    is a different and narrower claim than "KRAS-mutant lines need this gene".

    Parameters
    ----------
    driver_gene : str, optional
        Driver every line in the cohort must carry, e.g. "KRAS" (default).
    comutation_gene : str, optional
        Second gene whose status splits the cohort, e.g. "TP53" (default), "SMAD4", "CDKN2A".
    comutation_event : str, optional
        "mutation" (default), "deletion" or "amplification" for the second gene.
    driver_allele : str, optional
        Restrict the cohort to one driver allele, e.g. "G12D". None uses any driver mutation.
    cancer_type : str, optional
        Cancer context (default: "Pancreatic Cancer").
    data_lake_path, mutation_csv_path : str, optional
        DepMap directory and an optional genotype table overriding the default call source.
    top_n : int, optional
        Number of candidates printed (default: 20).
    p_threshold, fdr_threshold, min_effect_difference, max_mutant_mean_effect : optional
        Significance and effect-size gates.
    min_cohens_d : float, optional
        Require Cohen's d at or below this value (e.g. -0.8 for a large effect). None disables it.
    min_pct_mutant_dependent : float, optional
        Require this percentage of the altered arm to be dependent (gene effect < the depletion
        threshold). Separates a shift of the whole group from a shift driven by a few lines.
    max_pct_all_dependent : float, optional
        Reject genes depleted in more than this percentage of ALL screened lines, a stricter and
        more explicit control than the pan-essential flag alone.
    exclude_pan_essential : bool, optional
        Drop common-essential genes (default: True).
    output_csv_path : str, optional
        Write the full table to this CSV path.

    Returns
    -------
    str
        A research log with both arms, the ranked candidates, a check of whether each hit is specific
        to the co-altered group or general to the cancer type, QC warnings and a CANDIDATE_GENES line.

    """
    import numpy as np

    driver = driver_gene.strip().upper()
    partner = comutation_gene.strip().upper()
    mode = _STRATIFY_ALIASES.get(str(comutation_event).strip().lower()) if "_STRATIFY_ALIASES" in globals() else None
    mode = mode or str(comutation_event).strip().lower()
    if mode not in ("mutation", "deletion", "amplification"):
        return f"FAILURE: comutation_event='{comutation_event}' must be mutation, deletion or amplification."

    altered_word = {"mutation": "mutant", "deletion": "deleted", "amplification": "amplified"}[mode]
    control_word = {"mutation": "wild-type", "deletion": "intact", "amplification": "neutral"}[mode]
    log = [
        "=" * 78,
        f"CO-MUTATION STRATIFIED SL - {cancer_type}, all lines {driver}"
        + (f" {driver_allele}" if driver_allele else "") + f"-mutant, split by {partner} {mode}",
        f"Run at {datetime.now(tz=_UTC).strftime('%Y-%m-%d %H:%M UTC')}",
        "=" * 78,
        "",
    ]

    try:
        bundle = _load_depmap(data_lake_path)
    except FileNotFoundError as e:
        return f"FAILURE: {e}"
    gene_effect = bundle["gene_effect"]

    cohort, cohort_note = _select_cancer_models(bundle["model"], cancer_type)
    cohort = cohort[cohort["ModelID"].isin(gene_effect.index)]
    if len(cohort) == 0:
        return f"FAILURE: no {cancer_type} line with CRISPR data."

    try:
        driver_annotated, driver_source = _annotate_mutation_status(
            cohort, driver, bundle["data_lake_path"], mutation_csv_path
        )
    except (RuntimeError, ValueError) as e:
        return f"FAILURE: {e}"
    driver_mutant_ids, _, _ = _allele_groups(driver_annotated, driver_allele)
    driver_cohort = driver_annotated[driver_annotated["ModelID"].isin(driver_mutant_ids)]

    log.append("STEP 1 | Cohort - every line carries the driver")
    log.append(f"  {cohort_note}")
    log.append(f"  {cancer_type} lines with CRISPR data: {len(cohort)}")
    log.append(
        f"  {driver}{' ' + driver_allele if driver_allele else ''}-mutant: {len(driver_cohort)} "
        f"(this is the whole cohort from here on; there is no {driver} wild-type arm)"
    )

    try:
        if mode == "mutation":
            split, partner_source = _annotate_mutation_status(
                driver_cohort, partner, bundle["data_lake_path"], mutation_csv_path
            )
        else:
            split, partner_source = _annotate_copy_number_status(
                driver_cohort, partner, bundle["data_lake_path"], mode, mutation_csv_path
            )
    except (RuntimeError, ValueError) as e:
        return f"FAILURE: {e}"

    group_a = split.loc[split["MutationStatus"] == "MUT", "ModelID"].tolist()
    group_b = split.loc[split["MutationStatus"] == "WT", "ModelID"].tolist()
    unknown = split.loc[split["MutationStatus"] == "UNKNOWN", "ModelID"].tolist()

    log.append("")
    log.append(f"STEP 2 | Split by {partner} {mode}")
    log.append(f"  {partner} call source: {partner_source}")
    log.append(f"  {driver}-mutant + {partner}-{altered_word}: {len(group_a)} lines")
    log.append(f"  {driver}-mutant + {partner}-{control_word}: {len(group_b)} lines")
    log.append(f"  Not profiled for {partner} (excluded): {len(unknown)}")
    log.append(
        f"  A: {', '.join(split.loc[split['MutationStatus'] == 'MUT', 'StrippedCellLineName'].head(20))}"
    )
    log.append(
        f"  B: {', '.join(split.loc[split['MutationStatus'] == 'WT', 'StrippedCellLineName'].head(20))}"
    )

    if min(len(group_a), len(group_b)) < MIN_GROUP_SIZE:
        log.append("")
        log.append(
            f"FAILURE: one arm has fewer than {MIN_GROUP_SIZE} lines. Try another co-alteration - in PDAC, "
            "CDKN2A deletion and SMAD4 mutation give more balanced arms than TP53, whose wild-type arm is tiny."
        )
        return "\n".join(log)

    results = _welch_scan(gene_effect, group_a, group_b)
    selected = results[
        (results["p_one_sided"] < p_threshold)
        & (results["q_value"] < fdr_threshold)
        & (results["effect_difference"] <= min_effect_difference)
        & (results["mean_a"] <= max_mutant_mean_effect)
        & _apply_effect_filters(results, min_cohens_d, min_pct_mutant_dependent, max_pct_all_dependent)
    ]
    if exclude_pan_essential:
        selected = selected[~selected["pan_essential"]]
    selected = selected.sort_values(["effect_difference", "p_one_sided"]).head(top_n).reset_index(drop=True)

    log.append("")
    log.append(f"STEP 3 | Differential dependency ({partner}-{altered_word} vs {partner}-{control_word})")
    log.append(f"  Genes tested: {len(results)}")
    log.append(f"  p < {p_threshold}: {(results['p_one_sided'] < p_threshold).sum()} | "
               f"q < {fdr_threshold}: {(results['q_value'] < fdr_threshold).sum()}")
    log.append(f"  Surviving all gates: {len(selected)}")

    log.append("")
    log.append(f"TOP {len(selected)} CANDIDATES")
    log.append(
        f"{'gene':<12}{'A mean':>9}{'B mean':>9}{'delta':>8}{'p':>10}{'q':>9}"
        f"{'%A dep':>8}{'%B dep':>8}{'%all dep':>10}"
    )
    log.append("-" * 83)
    for _, row in selected.iterrows():
        log.append(
            f"{row['gene']:<12}{row['mean_a']:>9.3f}{row['mean_b']:>9.3f}{row['effect_difference']:>8.3f}"
            f"{row['p_one_sided']:>10.2e}{row['q_value']:>9.4f}{row['pct_a_dependent']:>7.0f}%"
            f"{row['pct_b_dependent']:>7.0f}%{row['pct_all_dependent']:>9.0f}%"
        )
    if len(selected) == 0:
        runners = results[
            (results["p_one_sided"] < p_threshold) & (results["effect_difference"] <= min_effect_difference)
        ]
        if exclude_pan_essential:
            runners = runners[~runners["pan_essential"]]
        log.append(f"  (none cleared the FDR gate; nominally significant runners-up, q >= {fdr_threshold}:)")
        for _, row in runners.nsmallest(min(top_n, 10), "effect_difference").iterrows():
            log.append(
                f"  {row['gene']:<12}A {row['mean_a']:>7.3f}  B {row['mean_b']:>7.3f}  "
                f"delta {row['effect_difference']:>7.3f}  p {row['p_one_sided']:.2e}  q {row['q_value']:.3f}"
            )

    log.append("")
    log.append(f"CANDIDATE_GENES: {', '.join(selected['gene'].tolist())}")
    log.append(
        f"  (these are dependencies of the {driver}-mutant / {partner}-{altered_word} genotype, not of "
        f"{driver} itself - a different and narrower claim)"
    )

    log.append("")
    log.append("QC WARNINGS")
    log.append(
        f"  - Both arms carry mutant {driver}, so nothing here is evidence about {driver} itself. The contrast "
        f"isolates what {partner} {mode} adds on top of it."
    )
    if min(len(group_a), len(group_b)) < INFORMATIVE_GROUP_SIZE:
        log.append(
            f"  - The smaller arm has {min(len(group_a), len(group_b))} lines (<{INFORMATIVE_GROUP_SIZE}); expect "
            "few or no genes to clear a genome-wide FDR. Shrink the gene set rather than raising the cutoff."
        )
    log.append(
        f"  - Co-alterations travel together: {partner} status correlates with other PDAC events, so a hit may "
        "track a co-occurring alteration rather than this one. check_dependency_confounders tests that."
    )
    if mode == "mutation":
        log.append(
            f"  - {partner} mutation calls do not separate loss-of-function from VUS. For TP53 in particular, a "
            "missense hotspot can be gain-of-function, which is a different biology from loss."
        )
    log.append(
        "  - The %all dep column is the whole DepMap panel. A gene depleted there as strongly as in arm A is a "
        "fitness gene, not a genotype-selective one."
    )

    if output_csv_path:
        os.makedirs(os.path.dirname(os.path.abspath(output_csv_path)), exist_ok=True)
        results.sort_values(["effect_difference", "p_one_sided"]).to_csv(output_csv_path, index=False)
        log.append(f"  Full table written to {output_csv_path}")

    log.append("")
    log.append("PROVENANCE")
    for item in bundle["provenance"]:
        log.append(f"  - {item}")
    log.append(f"  - {driver} calls: {driver_source}")
    log.append(f"  - {partner} calls: {partner_source}")
    log.append(
        f"  - Statistics: one-sided Welch t-test, BH-FDR across {len(results)} genes; depletion threshold "
        f"{DEPLETION_THRESHOLD}, pan-essential cutoff {PAN_ESSENTIAL_FRACTION:.0%}"
    )
    _ = np  # numpy is imported for the callers that extend this report
    return "\n".join(log)


# ---------------------------------------------------------------------------
# Allele-resolved comparison: G12D vs G12V
# ---------------------------------------------------------------------------
# Pathways where inhibiting KRAS G12D with MRTX1133 is expected to leave, or create, a second
# liability - the axes a combination arm would target. This is CURATED PRIOR KNOWLEDGE as of the
# date below, not something derived from the DepMap data, and it is used only to annotate
# candidates. A gene outside these sets is not thereby a worse candidate.
MRTX1133_COMBINATION_AS_OF = "2026-10"
MRTX1133_COMBINATION_PATHWAYS = {
    "RTK feedback reactivation": {
        "genes": {"EGFR", "ERBB2", "ERBB3", "FGFR1", "FGFR2", "IGF1R", "MET", "AXL", "GRB2", "SHC1"},
        "rationale": "KRAS inhibition relieves negative feedback and upstream RTKs re-activate RAS; "
        "RTK or pan-ERBB co-inhibition is the most advanced combination axis clinically.",
    },
    "Vertical RAS pathway": {
        "genes": {"SOS1", "PTPN11", "RAF1", "BRAF", "MAP2K1", "MAP2K2", "MAPK1", "MAPK3", "RASA1", "NF1", "SHOC2"},
        "rationale": "Hitting the same pathway above and below the driver suppresses the adaptive "
        "rebound that single-agent KRAS inhibition produces (SOS1 and SHP2 combinations are in trials).",
    },
    "PI3K-AKT-mTOR": {
        "genes": {"PIK3CA", "PIK3CB", "PIK3R1", "AKT1", "AKT2", "MTOR", "RICTOR", "RPTOR", "PDPK1"},
        "rationale": "The parallel effector arm of RAS; it sustains survival signalling when the MAPK arm is shut down.",
    },
    "Cell cycle re-entry": {
        "genes": {"CCND1", "CCND3", "CDK4", "CDK6", "CDK2", "CCNE1", "RB1", "SKP2"},
        "rationale": "KRAS inhibition arrests rather than kills; CDK4/6 co-inhibition deepens the arrest "
        "and is a tractable combination with launched drugs.",
    },
    "Apoptotic priming": {
        "genes": {"BCL2L1", "MCL1", "BCL2", "BAX", "BAK1", "BID", "PMAIP1", "BBC3"},
        "rationale": "A cytostatic response becomes cytotoxic when the mitochondrial apoptotic threshold "
        "is lowered; BCL-xL/MCL1 dependence is the usual bottleneck.",
    },
    "Autophagy and lysosome": {
        "genes": {"ATG7", "ATG5", "ATG3", "ULK1", "RAB7A", "LAMP1", "TFEB", "PIK3C3"},
        "rationale": "RAS-pathway inhibition induces protective autophagy in pancreatic cells; blocking it "
        "converts adaptation into death.",
    },
    "YAP/TAZ-TEAD bypass": {
        "genes": {"YAP1", "WWTR1", "TEAD1", "TEAD2", "TEAD3", "TEAD4", "AMOTL2", "LATS1", "LATS2"},
        "rationale": "The best-documented KRAS-independence route: tumours that escape KRAS inhibition "
        "re-route survival through YAP/TAZ-TEAD transcription.",
    },
}


def _combination_pathway(gene: str) -> tuple:
    """Return (pathway name, rationale) when ``gene`` sits on an MRTX1133 combination axis."""
    gene = gene.strip().upper()
    for name, spec in MRTX1133_COMBINATION_PATHWAYS.items():
        if gene in spec["genes"]:
            return name, spec["rationale"]
    return "", ""


def compare_allele_specific_dependencies(
    driver_gene: str = "KRAS",
    allele_a: str = "G12D",
    allele_b: str = "G12V",
    discovery_cohort: str = "pan-cancer",
    context_cancer_type: str = "Pancreatic Cancer",
    data_lake_path: str | None = None,
    mutation_csv_path: str | None = None,
    top_n: int = 15,
    p_threshold: float = 0.05,
    fdr_threshold: float = 0.25,
    min_effect_difference: float = -0.2,
    max_mutant_mean_effect: float = -0.3,
    exclude_pan_essential: bool = True,
    output_csv_path: str | None = None,
) -> str:
    """Separate dependencies shared by two driver alleles from those specific to one of them.

    "KRAS-mutant" is not one genotype. G12D and G12V differ in GTP hydrolysis rate, effector
    preference and - decisively for the clinic - in which inhibitors bind them, so a partner that
    only matters for one allele is lost when the alleles are pooled. This tool runs three contrasts
    rather than one:

    * allele A vs driver wild-type, and allele B vs driver wild-type - which partners each allele has
    * **allele A vs allele B directly** - the only contrast that can establish specificity, because a
      gene can clear the gate against wild-type in one allele and miss it in the other purely by power

    A candidate is reported A-SPECIFIC only when it passes against wild-type AND the direct A-vs-B
    contrast supports it; passing one but not the other is reported as A-ONLY (unconfirmed), which is
    a weaker claim. Candidates are then checked for presence in ``context_cancer_type`` lines, and
    annotated with the MRTX1133 combination axis they sit on, where they sit on one.

    Parameters
    ----------
    driver_gene : str, optional
        Driver carrying the alleles (default: "KRAS").
    allele_a, allele_b : str, optional
        Protein changes to compare, e.g. "G12D" (default) and "G12V".
    discovery_cohort : str, optional
        "pan-cancer" (default) or a cancer type. Pan-cancer is usually required: within PDAC the
        driver wild-type arm is too small for either allele contrast.
    context_cancer_type : str, optional
        Cancer type candidates are checked for presence in (default: "Pancreatic Cancer").
    data_lake_path, mutation_csv_path : str, optional
        DepMap directory and an optional genotype table overriding the default call source.
    top_n : int, optional
        Number of rows printed per class (default: 15).
    p_threshold, fdr_threshold, min_effect_difference, max_mutant_mean_effect : optional
        Significance and effect-size gates applied to each allele-vs-wild-type scan.
    exclude_pan_essential : bool, optional
        Drop common-essential genes (default: True).
    output_csv_path : str, optional
        Write the joined three-contrast table to this CSV path.

    Returns
    -------
    str
        A research log with the three cohorts, allele-A-specific candidates first (clinical priority),
        then allele-B-specific, then shared, each with the direct contrast, the context check and the
        combination-pathway annotation.

    """
    import numpy as np
    import pandas as pd

    driver = driver_gene.strip().upper()
    allele_a, allele_b = allele_a.strip().upper(), allele_b.strip().upper()
    log = [
        "=" * 78,
        f"ALLELE-RESOLVED DEPENDENCIES - {driver} {allele_a} vs {allele_b}",
        f"Run at {datetime.now(tz=_UTC).strftime('%Y-%m-%d %H:%M UTC')}",
        "=" * 78,
        "",
    ]

    try:
        bundle = _load_depmap(data_lake_path)
    except FileNotFoundError as e:
        return f"FAILURE: {e}"
    gene_effect = bundle["gene_effect"]

    cohort, cohort_note = _select_cancer_models(bundle["model"], discovery_cohort)
    cohort = cohort[cohort["ModelID"].isin(gene_effect.index)]
    try:
        annotated, source = _annotate_mutation_status(cohort, driver, bundle["data_lake_path"], mutation_csv_path)
    except (RuntimeError, ValueError) as e:
        return f"FAILURE: {e}"

    ids_a, wildtype_ids, _ = _allele_groups(annotated, allele_a)
    ids_b, _, _ = _allele_groups(annotated, allele_b)
    if min(len(ids_a), len(ids_b)) < MIN_GROUP_SIZE or len(wildtype_ids) < MIN_GROUP_SIZE:
        return (
            f"FAILURE: group sizes {allele_a}={len(ids_a)}, {allele_b}={len(ids_b)}, "
            f"wild-type={len(wildtype_ids)}. Use discovery_cohort='pan-cancer' or check the allele spelling."
        )

    log.append("STEP 1 | Cohorts")
    log.append(f"  {cohort_note}")
    log.append(f"  Genotype source: {source}")
    log.append(f"  {driver} {allele_a}: {len(ids_a)} lines | {driver} {allele_b}: {len(ids_b)} lines "
               f"| {driver} wild-type: {len(wildtype_ids)} lines")
    lineages_a = annotated[annotated["ModelID"].isin(ids_a)]["OncotreeLineage"].value_counts()
    lineages_b = annotated[annotated["ModelID"].isin(ids_b)]["OncotreeLineage"].value_counts()
    log.append(f"  {allele_a} lineages: {lineages_a.head(4).to_dict()}")
    log.append(f"  {allele_b} lineages: {lineages_b.head(4).to_dict()}")

    scan_a = _welch_scan(gene_effect, ids_a, wildtype_ids).set_index("gene")
    scan_b = _welch_scan(gene_effect, ids_b, wildtype_ids).set_index("gene")
    direct = _welch_scan(gene_effect, ids_a, ids_b).set_index("gene")

    def passes(frame):
        ok = (
            (frame["p_one_sided"] < p_threshold)
            & (frame["q_value"] < fdr_threshold)
            & (frame["effect_difference"] <= min_effect_difference)
            & (frame["mean_a"] <= max_mutant_mean_effect)
        )
        if exclude_pan_essential:
            ok &= ~frame["pan_essential"]
        return ok

    pass_a, pass_b = passes(scan_a), passes(scan_b)
    log.append("")
    log.append(f"STEP 2 | Two scans against the same {driver} wild-type arm")
    log.append(f"  {allele_a} vs wild-type: {int(pass_a.sum())} genes clear every gate")
    log.append(f"  {allele_b} vs wild-type: {int(pass_b.sum())} genes clear every gate")

    genes = sorted(set(scan_a.index[pass_a]) | set(scan_b.index[pass_b]))
    if not genes:
        log.append("")
        log.append("No gene cleared the gates for either allele. Loosen fdr_threshold or restrict the gene set.")
        return "\n".join(log)

    # Context presence: are these dependencies actually there in the cancer type of interest?
    context_cohort, _ = _select_cancer_models(bundle["model"], context_cancer_type)
    context_cohort = context_cohort[context_cohort["ModelID"].isin(gene_effect.index)]
    context_annotated, _ = _annotate_mutation_status(
        context_cohort, driver, bundle["data_lake_path"], mutation_csv_path
    )
    context_a, _, _ = _allele_groups(context_annotated, allele_a)
    context_b, _, _ = _allele_groups(context_annotated, allele_b)
    log.append("")
    log.append(f"STEP 3 | Context: {context_cancer_type} - {allele_a} {len(context_a)} lines, "
               f"{allele_b} {len(context_b)} lines")

    rows = []
    for gene in genes:
        in_a, in_b = gene in scan_a.index, gene in scan_b.index
        a_pass = bool(pass_a.get(gene, False))
        b_pass = bool(pass_b.get(gene, False))
        direct_row = direct.loc[gene] if gene in direct.index else None
        direct_delta = float(direct_row["effect_difference"]) if direct_row is not None else np.nan
        direct_p = float(direct_row["p_one_sided"]) if direct_row is not None else np.nan
        # The direct contrast is one-sided for "A more depleted"; flip it to read B-specificity.
        direct_p_b = 1.0 - direct_p if np.isfinite(direct_p) else np.nan
        confirms_a = np.isfinite(direct_delta) and direct_delta <= min_effect_difference and direct_p < p_threshold
        confirms_b = np.isfinite(direct_delta) and direct_delta >= -min_effect_difference and direct_p_b < p_threshold

        if a_pass and b_pass:
            klass = "SHARED"
        elif a_pass and confirms_a:
            klass = f"{allele_a}-SPECIFIC"
        elif a_pass:
            klass = f"{allele_a}-ONLY (direct contrast does not confirm)"
        elif b_pass and confirms_b:
            klass = f"{allele_b}-SPECIFIC"
        else:
            klass = f"{allele_b}-ONLY (direct contrast does not confirm)"

        context_ids = context_a if klass.startswith(allele_a) or klass == "SHARED" else context_b
        column = bundle["gene_columns"].get(gene)
        context_mean, context_pct = np.nan, np.nan
        if column is not None and len(context_ids) >= MIN_GROUP_SIZE:
            values = gene_effect.loc[context_ids, column].dropna()
            if len(values):
                context_mean = float(values.mean())
                context_pct = float((values < DEPLETION_THRESHOLD).mean() * 100)
        pathway, rationale = _combination_pathway(gene)
        rows.append(
            {
                "gene": gene,
                "class": klass,
                "delta_a_vs_wt": float(scan_a.loc[gene, "effect_difference"]) if in_a else np.nan,
                "q_a_vs_wt": float(scan_a.loc[gene, "q_value"]) if in_a else np.nan,
                "delta_b_vs_wt": float(scan_b.loc[gene, "effect_difference"]) if in_b else np.nan,
                "q_b_vs_wt": float(scan_b.loc[gene, "q_value"]) if in_b else np.nan,
                "delta_a_vs_b": direct_delta,
                "p_a_vs_b": direct_p,
                "context_mean_effect": context_mean,
                "context_pct_dependent": context_pct,
                "combination_pathway": pathway,
                "combination_rationale": rationale,
                "pct_all_dependent": float(scan_a.loc[gene, "pct_all_dependent"]) if in_a
                else float(scan_b.loc[gene, "pct_all_dependent"]),
            }
        )
    table = pd.DataFrame(rows)

    def block(title, subset, note):
        log.append("")
        log.append(f"{title}  ({len(subset)})")
        if len(subset) == 0:
            log.append("  none")
            return
        log.append(f"  {note}")
        log.append(
            f"  {'gene':<10}{'dA/WT':>8}{'qA':>8}{'dB/WT':>8}{'qB':>8}{'dA-B':>8}{'pA-B':>9}"
            f"{'ctx':>8}{'ctx%':>7}  combination axis"
        )
        log.append("  " + "-" * 88)
        for _, row in subset.head(top_n).iterrows():
            fmt = lambda v, d=3: "  n/a" if not np.isfinite(v) else f"{v:.{d}f}"
            context_pct = "n/a" if not np.isfinite(row["context_pct_dependent"]) else f"{row['context_pct_dependent']:.0f}%"
            log.append(
                f"  {row['gene']:<10}{fmt(row['delta_a_vs_wt']):>8}{fmt(row['q_a_vs_wt'], 4):>8}"
                f"{fmt(row['delta_b_vs_wt']):>8}{fmt(row['q_b_vs_wt'], 4):>8}{fmt(row['delta_a_vs_b']):>8}"
                f"{fmt(row['p_a_vs_b'], 4):>9}{fmt(row['context_mean_effect']):>8}{context_pct:>7}"
                f"  {row['combination_pathway'] or '-'}"
            )

    specific_a = table[table["class"] == f"{allele_a}-SPECIFIC"].sort_values("delta_a_vs_b")
    only_a = table[table["class"].str.startswith(f"{allele_a}-ONLY")].sort_values("delta_a_vs_wt")
    specific_b = table[table["class"] == f"{allele_b}-SPECIFIC"].sort_values("delta_a_vs_b", ascending=False)
    only_b = table[table["class"].str.startswith(f"{allele_b}-ONLY")].sort_values("delta_b_vs_wt")
    shared = table[table["class"] == "SHARED"].sort_values("delta_a_vs_wt")

    log.append("")
    log.append("=" * 78)
    log.append(f"CLINICAL PRIORITY 1 - {allele_a}-SPECIFIC (direct {allele_a} vs {allele_b} contrast confirms)")
    log.append("=" * 78)
    block(f"{allele_a}-SPECIFIC", specific_a,
          f"stronger in {allele_a} than in both wild-type AND {allele_b} - the only class where allele "
          "specificity is actually demonstrated rather than inferred from two separate tests")
    block(f"{allele_a}-ONLY (unconfirmed)", only_a,
          f"clears the gate against wild-type for {allele_a} but the direct contrast against {allele_b} does "
          "not separate them; may be a power difference between the two arms, not biology")
    block("SHARED by both alleles", shared,
          f"a partner of {driver} mutation generally; a broader patient population but no allele rationale")
    block(f"{allele_b}-SPECIFIC", specific_b, f"stronger in {allele_b} - listed for completeness")
    block(f"{allele_b}-ONLY (unconfirmed)", only_b, "")

    combo = table[(table["combination_pathway"] != "") & table["class"].str.startswith(allele_a)]
    log.append("")
    log.append(f"MRTX1133 COMBINATION AXES among {allele_a} candidates  ({len(combo)})")
    if len(combo) == 0:
        log.append(
            f"  No {allele_a} candidate falls on a curated MRTX1133 combination axis. That is a statement about "
            "the curated pathway sets, not evidence against combination potential."
        )
    for pathway in sorted(combo["combination_pathway"].unique()):
        members = combo[combo["combination_pathway"] == pathway]
        log.append(f"  [{pathway}] {', '.join(members['gene'])}")
        log.append(f"      {MRTX1133_COMBINATION_PATHWAYS[pathway]['rationale']}")

    log.append("")
    log.append(f"CANDIDATE_GENES: {', '.join(specific_a['gene'].head(top_n))}")
    log.append(f"ALLELE_SPECIFIC_{allele_a}: {', '.join(specific_a['gene'])}")
    log.append(f"ALLELE_SPECIFIC_{allele_b}: {', '.join(specific_b['gene'])}")
    log.append(f"SHARED_PARTNERS: {', '.join(shared['gene'].head(top_n))}")

    log.append("")
    log.append("QC WARNINGS")
    log.append(
        f"  - The two alleles have different arm sizes ({len(ids_a)} vs {len(ids_b)}), so the larger arm clears "
        "gates more easily. That is why a gene passing for one allele only is reported as unconfirmed unless the "
        "direct contrast separates them."
    )
    log.append(
        f"  - Lineage composition differs between the arms ({allele_a}: {lineages_a.index[0]} "
        f"{lineages_a.iloc[0] / len(ids_a):.0%}; {allele_b}: {lineages_b.index[0]} "
        f"{lineages_b.iloc[0] / len(ids_b):.0%}). An allele-specific hit can be a lineage-specific hit."
    )
    log.append(
        f"  - Combination axes are curated prior knowledge as of {MRTX1133_COMBINATION_AS_OF}, not a result from "
        "this data. A gene outside those sets is not a worse candidate; it simply has no annotated rationale here."
    )
    log.append(
        "  - Allele identity comes from the variant strings of the mutation source; a line with two KRAS variants "
        "is counted for each allele it carries."
    )

    if output_csv_path:
        os.makedirs(os.path.dirname(os.path.abspath(output_csv_path)), exist_ok=True)
        table.to_csv(output_csv_path, index=False)
        log.append(f"  Full table written to {output_csv_path} ({len(table)} genes)")

    log.append("")
    log.append("PROVENANCE")
    for item in bundle["provenance"]:
        log.append(f"  - {item}")
    log.append(f"  - Genotype calls: {source}")
    log.append(
        f"  - Statistics: three one-sided Welch scans ({allele_a} vs WT, {allele_b} vs WT, {allele_a} vs "
        f"{allele_b}), BH-FDR within each scan; depletion threshold {DEPLETION_THRESHOLD}"
    )
    return "\n".join(log)


# ---------------------------------------------------------------------------
# Biological validation of candidates
# ---------------------------------------------------------------------------
STRING_API = "https://string-db.org/api"
UNIPROT_API = "https://rest.uniprot.org/uniprotkb/search"
ALPHAFOLD_API = "https://alphafold.ebi.ac.uk/api/prediction"
TCGA_PAAD_STUDY = "paad_tcga_pan_can_atlas_2018"

# Which KRAS effector axis a gene sits on. CURATED PRIOR KNOWLEDGE - an axis label is a hypothesis
# about mechanism, not a measurement from this screen.
KRAS_AXIS_AS_OF = "2026-10"
KRAS_EFFECTOR_AXES = {
    "MAPK/ERK": {"RAF1", "BRAF", "ARAF", "MAP2K1", "MAP2K2", "MAPK1", "MAPK3", "SHOC2", "KSR1", "DUSP4",
                 "SPRY2", "ETV4", "ETV5", "KLF5", "MYC", "FOSL1", "CCND1", "SKP2"},
    "PI3K": {"PIK3CA", "PIK3CB", "PIK3R1", "AKT1", "AKT2", "MTOR", "RICTOR", "RPTOR", "PDPK1", "PTEN", "RHEB"},
    "RAL": {"RALA", "RALB", "RALGDS", "RGL1", "RGL2", "EXOC1", "EXOC2", "EXOC4", "TBK1"},
    "RAC1/cytoskeleton": {"RAC1", "DOCK5", "DOCK1", "DOCK2", "TIAM1", "PAK1", "PAK4", "ARHGEF7", "CDC42", "WASF2"},
    "Macropinocytosis/trafficking": {"RAB10", "RAB7A", "RAB5A", "PIKFYVE", "NHE1", "SLC9A1", "ATG7", "LAMP1",
                                     "SNAP23", "VAMP7", "EXOC3", "AP2M1", "CTNS", "SLC38A9"},
    "Metabolic rewiring": {"SLC7A5", "SLC7A11", "GOT1", "GLS", "SREBF1", "SCAP", "ACLY", "FASN", "CHKA", "PHGDH"},
    "WNT/transcription": {"CTNNB1", "TCF7L2", "LEF1", "TCF7", "YAP1", "WWTR1", "TEAD1"},
}


def _kras_axis(gene: str) -> str:
    gene = gene.strip().upper()
    for axis, members in KRAS_EFFECTOR_AXES.items():
        if gene in members:
            return axis
    return ""


def _pdac_restricted_effect(gene: str, bundle: dict, allele: str, mutation_csv_path: str | None) -> dict:
    """Cohen's d for the gene inside PDAC lines only, against two different control arms.

    The protocol asks for a PDAC-restricted effect size, but PDAC has only a handful of KRAS
    wild-type lines, so that contrast alone is not interpretable. The better-powered contrast -
    allele lines against every other PDAC line - is reported beside it.
    """
    import numpy as np

    column = bundle["gene_columns"].get(gene.strip().upper())
    if column is None:
        return {}
    cohort, _ = _select_cancer_models(bundle["model"], "Pancreatic Cancer")
    cohort = cohort[cohort["ModelID"].isin(bundle["gene_effect"].index)]
    annotated, _ = _annotate_mutation_status(cohort, "KRAS", bundle["data_lake_path"], mutation_csv_path)
    allele_ids, wildtype_ids, other_ids = _allele_groups(annotated, allele)

    def effect(ids_a, ids_b):
        a = bundle["gene_effect"].loc[ids_a, column].dropna()
        b = bundle["gene_effect"].loc[ids_b, column].dropna()
        if len(a) < 2 or len(b) < 2:
            return {"n_a": len(a), "n_b": len(b), "cohens_d": np.nan, "delta": np.nan}
        pooled = np.sqrt((a.std(ddof=1) ** 2 + b.std(ddof=1) ** 2) / 2)
        return {
            "n_a": len(a), "n_b": len(b), "delta": float(a.mean() - b.mean()),
            "cohens_d": float((a.mean() - b.mean()) / pooled) if pooled > 0 else np.nan,
        }

    return {"vs_wildtype": effect(allele_ids, wildtype_ids), "vs_other_pdac": effect(allele_ids, wildtype_ids + other_ids)}


def _logrank(times_a, events_a, times_b, events_b) -> tuple:
    """Two-sample log-rank test. Returns (chi2, p, observed_a, expected_a)."""
    import numpy as np
    from scipy import stats

    times = np.concatenate([times_a, times_b])
    events = np.concatenate([events_a, events_b])
    group = np.concatenate([np.ones(len(times_a)), np.zeros(len(times_b))])
    observed_a = float(events[group == 1].sum())

    expected_a, variance = 0.0, 0.0
    for time in np.unique(times[events == 1]):
        at_risk = times >= time
        n_total, n_a = at_risk.sum(), (at_risk & (group == 1)).sum()
        d_total = ((times == time) & (events == 1)).sum()
        if n_total <= 1:
            continue
        expected_a += d_total * n_a / n_total
        variance += d_total * (n_total - d_total) * n_a * (n_total - n_a) / (n_total**2 * (n_total - 1))
    if variance <= 0:
        return np.nan, np.nan, observed_a, expected_a
    chi2 = (observed_a - expected_a) ** 2 / variance
    return float(chi2), float(stats.chi2.sf(chi2, 1)), observed_a, float(expected_a)


def _tcga_paad_survival(gene: str) -> dict:
    """Split TCGA PAAD patients at the median expression of ``gene`` and log-rank their survival."""
    import numpy as np

    clinical = {}
    for attribute in ("OS_MONTHS", "OS_STATUS"):
        payload = _http_get(
            f"{CBIOPORTAL_API}/studies/{TCGA_PAAD_STUDY}/clinical-data",
            params={"attributeId": attribute, "clinicalDataType": "PATIENT", "pageSize": 1000},
            timeout=60, retries=2,
        )
        if payload is None:
            return {"status": "unavailable", "note": "cBioPortal clinical data could not be retrieved"}
        for record in payload:
            clinical.setdefault(record["patientId"], {})[attribute] = record["value"]

    gene_info = _http_get(f"{CBIOPORTAL_API}/genes/{urllib.parse.quote(gene.upper())}", timeout=30, retries=2)
    if not gene_info or "entrezGeneId" not in gene_info:
        return {"status": "unavailable", "note": f"{gene} not resolvable in cBioPortal"}

    expression = _http_post(
        f"{CBIOPORTAL_API}/molecular-profiles/{TCGA_PAAD_STUDY}_rna_seq_v2_mrna_median_Zscores/molecular-data/fetch",
        payload={"entrezGeneIds": [gene_info["entrezGeneId"]], "sampleListId": f"{TCGA_PAAD_STUDY}_rna_seq_v2_mrna"},
        params={"projection": "SUMMARY"}, timeout=60, retries=2,
    )
    if not expression:
        return {"status": "unavailable", "note": "expression profile could not be retrieved"}

    rows = []
    for record in expression:
        patient = record.get("patientId") or str(record.get("sampleId", ""))[:12]
        info = clinical.get(patient)
        value = record.get("value")
        if not info or value is None or "OS_MONTHS" not in info or "OS_STATUS" not in info:
            continue
        try:
            months = float(info["OS_MONTHS"])
        except (TypeError, ValueError):
            continue
        rows.append((float(value), months, 1 if info["OS_STATUS"].startswith("1") else 0))
    if len(rows) < 30:
        return {"status": "unavailable", "note": f"only {len(rows)} patients with both expression and survival"}

    values = np.array([r[0] for r in rows])
    times = np.array([r[1] for r in rows])
    events = np.array([r[2] for r in rows])
    high = values > np.median(values)
    chi2, p, observed, expected = _logrank(times[high], events[high], times[~high], events[~high])
    median_high = float(np.median(times[high]))
    median_low = float(np.median(times[~high]))
    if not np.isfinite(p):
        return {"status": "unavailable", "note": "log-rank variance was zero"}
    worse = observed > expected
    return {
        "status": "computed",
        "n_patients": len(rows),
        "n_high": int(high.sum()),
        "median_os_high": median_high,
        "median_os_low": median_low,
        "logrank_p": p,
        "verdict": ("high expression → worse OS" if worse else "high expression → better OS")
        + (" (significant)" if p < 0.05 else " (not significant)"),
        "answer": "yes" if (worse and p < 0.05) else ("no" if p < 0.05 else "not significant"),
    }


def _string_neighbours(gene: str, min_score: float = 0.7, limit: int = 20) -> list:
    """First-shell STRING partners of ``gene`` above a combined-score cutoff."""
    payload = _http_get(
        f"{STRING_API}/json/network",
        params={"identifiers": gene, "species": 9606, "required_score": int(min_score * 1000)},
        timeout=45, retries=2,
    )
    if not payload:
        return []
    neighbours = {}
    for edge in payload:
        for side, other in (("preferredName_A", "preferredName_B"), ("preferredName_B", "preferredName_A")):
            if str(edge.get(side, "")).upper() == gene.upper():
                partner = str(edge.get(other, "")).upper()
                if partner and partner != gene.upper():
                    neighbours[partner] = max(neighbours.get(partner, 0), float(edge.get("score", 0)))
    return sorted(neighbours.items(), key=lambda kv: -kv[1])[:limit]


def _chembl_target_class(gene: str) -> dict:
    """ChEMBL target type and protein class for ``gene``; an empty dict when unavailable."""
    payload = _http_get(
        f"{CHEMBL_API}/target.json",
        params={"target_components__target_component_synonyms__component_synonym": gene, "limit": 10},
        timeout=20, retries=1,
    )
    if payload is None:
        return {"status": "unavailable"}
    for target in payload.get("targets", []):
        if target.get("organism") == "Homo sapiens" and target.get("target_type") == "SINGLE PROTEIN":
            classes = []
            for component in target.get("target_components", []):
                for entry in component.get("target_component_xrefs", []):
                    if entry.get("xref_src_db") == "ChEMBL_Protein_Class":
                        classes.append(entry.get("xref_name", ""))
            return {"status": "found", "target_chembl_id": target["target_chembl_id"],
                    "target_type": target["target_type"], "protein_class": "; ".join(filter(None, classes))}
    return {"status": "not found"}


def _alphafold_confidence(gene: str) -> dict:
    """Mean pLDDT of the AlphaFold model for ``gene``.

    This is model confidence, NOT pocket quality: a pocket-quality score needs a pocket detector
    (fpocket, PocketMiner) run over the structure, which this environment does not have. Low pLDDT
    does rule the gene out of structure-based pocket work, so the number is still decision-relevant.
    """
    import numpy as np

    search = _http_get(
        UNIPROT_API,
        params={"query": f"gene_exact:{gene} AND organism_id:9606 AND reviewed:true",
                "fields": "accession", "format": "json", "size": 1},
        timeout=30, retries=2,
    )
    accessions = [r["primaryAccession"] for r in (search or {}).get("results", [])]
    if not accessions:
        return {"status": "unavailable", "note": "no reviewed UniProt entry"}

    prediction = _http_get(f"{ALPHAFOLD_API}/{accessions[0]}", timeout=45, retries=2)
    if not prediction:
        return {"status": "unavailable", "note": "AlphaFold DB has no model", "uniprot": accessions[0]}
    url = prediction[0].get("pdbUrl")
    text = _http_get(url, timeout=60, retries=1, as_json=False) if url else None
    if not text:
        return {"status": "unavailable", "note": "model file not retrievable", "uniprot": accessions[0]}

    scores = [float(line[60:66]) for line in text.splitlines()
              if line.startswith("ATOM") and line[12:16].strip() == "CA"]
    if not scores:
        return {"status": "unavailable", "note": "no CA atoms parsed", "uniprot": accessions[0]}
    mean = float(np.mean(scores))
    band = "very high (>90)" if mean > 90 else "confident (70-90)" if mean > 70 else "low (<70) - likely disordered"
    return {"status": "computed", "uniprot": accessions[0], "mean_plddt": round(mean, 1),
            "n_residues": len(scores), "band": band}


def validate_candidates_biologically(
    candidate_genes,
    driver_gene: str = "KRAS",
    allele: str = "G12D",
    data_lake_path: str | None = None,
    mutation_csv_path: str | None = None,
    string_min_score: float = 0.7,
    max_neighbours: int = 15,
    include_survival: bool = True,
    include_structure: bool = True,
    output_csv_path: str | None = None,
) -> str:
    """Run the biological validation battery on candidate genes, one evidence axis at a time.

    Statistical selection says a dependency is real in cell lines. None of these checks repeat that;
    each asks a different question that cell-line statistics cannot answer:

    * **PDAC-restricted effect** - does the effect survive inside the cancer type, where the control
      arm is small? Reported against both the wild-type arm and every other PDAC line.
    * **KRAS effector axis** - which downstream arm the gene sits on, as a mechanistic hypothesis.
    * **TCGA PAAD survival** - does expression track outcome in patients? A log-rank test on the
      median split, which is clinical association, not causality.
    * **Drug status** - a direct inhibitor, or a druggable STRING neighbour within one hop, which is
      how an undruggable node still becomes an experiment.
    * **Target class and structure confidence** - ChEMBL target class, plus mean AlphaFold pLDDT.
      pLDDT is model confidence and a prerequisite for pocket work; it is NOT a pocket-quality score,
      which needs a pocket detector this environment does not have.

    Parameters
    ----------
    candidate_genes : list[str] | str
        Genes to validate.
    driver_gene, allele : str, optional
        Genotype defining the PDAC-restricted contrast (default: KRAS G12D).
    data_lake_path, mutation_csv_path : str, optional
        DepMap directory and an optional genotype table.
    string_min_score : float, optional
        STRING combined-score cutoff for the one-hop neighbourhood (default: 0.7).
    max_neighbours : int, optional
        Neighbours examined per gene when looking for a druggable one (default: 15).
    include_survival : bool, optional
        Run the TCGA PAAD log-rank test (default: True; needs cBioPortal).
    include_structure : bool, optional
        Fetch AlphaFold model confidence (default: True; needs UniProt and AlphaFold DB).
    output_csv_path : str, optional
        Write the per-candidate table to this CSV path.

    Returns
    -------
    str
        A per-candidate validation report plus a summary table, and a separate list of druggable
        network neighbours for candidates that have no direct inhibitor.

    """
    import numpy as np
    import pandas as pd

    genes = _parse_gene_list(candidate_genes)
    if not genes:
        return "FAILURE: no candidate gene was supplied."
    try:
        bundle = _load_depmap(data_lake_path)
    except FileNotFoundError as e:
        return f"FAILURE: {e}"

    hub = _load_repurposing_hub(data_lake_path)
    hub_table = hub["table"] if hub is not None else None

    def inhibitors_for(gene):
        drugs = _drugs_from_curated(gene)
        if hub_table is not None:
            drugs += hub_table[hub_table["gene"] == gene].to_dict("records")
        found = _drugs_from_chembl(gene)
        if found:
            drugs += found
        return [d for d in drugs if _is_inhibitory(d.get("moa", ""))]

    log = [
        "=" * 78,
        f"BIOLOGICAL VALIDATION - {len(genes)} candidate(s), {driver_gene.upper()} {allele} context",
        f"Run at {datetime.now(tz=_UTC).strftime('%Y-%m-%d %H:%M UTC')}",
        "=" * 78,
    ]
    rows, neighbour_options = [], {}
    for gene in genes:
        log.append("")
        log.append(f"### {gene}")
        effect = _pdac_restricted_effect(gene, bundle, allele, mutation_csv_path)
        if effect:
            wt, other = effect["vs_wildtype"], effect["vs_other_pdac"]
            log.append(
                f"  PDAC-restricted effect: vs {driver_gene.upper()} wild-type d={wt['cohens_d']:.2f} "
                f"(n={wt['n_a']}v{wt['n_b']}) | vs all other PDAC lines d={other['cohens_d']:.2f} "
                f"(n={other['n_a']}v{other['n_b']})"
            )
        axis = _kras_axis(gene)
        log.append(f"  KRAS effector axis: {axis or 'not assigned in the curated map'}")

        survival = _tcga_paad_survival(gene) if include_survival else {"status": "skipped"}
        if survival.get("status") == "computed":
            log.append(
                f"  TCGA PAAD survival: {survival['verdict']} | median OS {survival['median_os_high']:.1f} vs "
                f"{survival['median_os_low']:.1f} months, log-rank p={survival['logrank_p']:.2e} "
                f"(n={survival['n_patients']})"
            )
        else:
            log.append(f"  TCGA PAAD survival: not available ({survival.get('note', survival.get('status'))})")

        direct = inhibitors_for(gene)
        clinical = [d for d in direct if d["phase_rank"] >= 2]
        if direct:
            best = max(direct, key=lambda d: d["phase_rank"])
            log.append(f"  Direct inhibitor: {len(direct)} ({len(clinical)} clinical) - "
                       f"{best['drug']} [{best['clinical_phase']}]")
            neighbours = []
        else:
            log.append("  Direct inhibitor: none")
            neighbours = []
            for partner, score in _string_neighbours(gene, string_min_score, max_neighbours):
                partner_drugs = [d for d in inhibitors_for(partner) if d["phase_rank"] >= 2]
                if partner_drugs:
                    best = max(partner_drugs, key=lambda d: d["phase_rank"])
                    neighbours.append({"neighbour": partner, "string_score": round(score, 3),
                                       "drug": best["drug"], "phase": best["clinical_phase"]})
            neighbour_options[gene] = neighbours
            if neighbours:
                log.append(f"  Druggable STRING neighbours (score >= {string_min_score}): " +
                           ", ".join(f"{n['neighbour']} ({n['drug']})" for n in neighbours[:4]))
            else:
                log.append(f"  Druggable STRING neighbours (score >= {string_min_score}): none found")

        target_class = _chembl_target_class(gene)
        log.append(f"  ChEMBL target class: {target_class.get('protein_class') or target_class.get('status')}")
        structure = _alphafold_confidence(gene) if include_structure else {"status": "skipped"}
        if structure.get("status") == "computed":
            log.append(f"  AlphaFold model confidence: mean pLDDT {structure['mean_plddt']} - {structure['band']} "
                       f"({structure['n_residues']} residues, {structure['uniprot']}). NOT a pocket score.")
        else:
            log.append(f"  AlphaFold model confidence: not available ({structure.get('note', structure.get('status'))})")

        rows.append({
            "gene": gene,
            "pdac_d_vs_wt": effect.get("vs_wildtype", {}).get("cohens_d", np.nan) if effect else np.nan,
            "pdac_d_vs_other": effect.get("vs_other_pdac", {}).get("cohens_d", np.nan) if effect else np.nan,
            "kras_axis": axis,
            "tcga_survival": survival.get("answer", "not available"),
            "tcga_logrank_p": survival.get("logrank_p", np.nan),
            "n_direct_inhibitors": len(direct),
            "n_clinical_inhibitors": len(clinical),
            "druggable_neighbours": "; ".join(f"{n['neighbour']}:{n['drug']}" for n in neighbours),
            "chembl_class": target_class.get("protein_class", ""),
            "mean_plddt": structure.get("mean_plddt", np.nan),
        })

    table = pd.DataFrame(rows)
    log.append("")
    log.append("SUMMARY")
    log.append(f"{'gene':<10}{'PDAC d(WT)':>11}{'PDAC d(oth)':>12}{'axis':<30}{'TCGA OS':>14}{'drug':>8}")
    log.append("-" * 86)
    for _, row in table.iterrows():
        log.append(
            f"{row['gene']:<10}{row['pdac_d_vs_wt']:>11.2f}{row['pdac_d_vs_other']:>12.2f}"
            f"{(row['kras_axis'] or '-'):<30}{row['tcga_survival']:>14}"
            f"{(str(row['n_clinical_inhibitors']) + ' clin' if row['n_clinical_inhibitors'] else 'none'):>8}"
        )

    no_drug = [g for g, n in neighbour_options.items() if n]
    if no_drug:
        log.append("")
        log.append("NO-DRUG CANDIDATES - druggable network neighbours within one hop")
        for gene in no_drug:
            log.append(f"  {gene}:")
            for option in neighbour_options[gene][:6]:
                log.append(f"    {option['neighbour']:<10} STRING {option['string_score']:.3f}  "
                           f"{option['drug']} [{option['phase']}]")

    log.append("")
    log.append("QC WARNINGS")
    log.append(
        "  - TCGA survival is an association in bulk tumour expression; it mixes tumour and stroma and says "
        "nothing about whether the gene is required. A gene can predict outcome without being a target."
    )
    log.append(
        "  - A druggable neighbour is a hypothesis about the pathway, not about the gene: inhibiting the "
        "neighbour tests a different node, and the knockout remains the direct test."
    )
    log.append(
        f"  - KRAS effector axis assignments are curated as of {KRAS_AXIS_AS_OF}, not derived from this screen."
    )
    log.append(
        "  - Mean pLDDT is AlphaFold's own confidence, NOT pocket quality. A pocket-quality score requires a "
        "pocket detector (fpocket, PocketMiner) over the structure, which is not available here."
    )

    if output_csv_path:
        os.makedirs(os.path.dirname(os.path.abspath(output_csv_path)), exist_ok=True)
        table.to_csv(output_csv_path, index=False)
        log.append(f"  Table written to {output_csv_path}")
    return "\n".join(log)
