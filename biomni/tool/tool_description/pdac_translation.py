description = [
    {
        "description": "Report which pancreatic cancer (PDAC) driver alterations actually split the DepMap cell "
        "line cohort into two testable groups, before a synthetic lethality run is committed to a driver. Counts "
        "mutant/wild-type (or amplified/neutral, deleted/neutral) lines for each curated PDAC driver and labels "
        "the contrast TESTABLE, UNDERPOWERED or UNTESTABLE - KRAS is altered in most PDAC lines, so its "
        "wild-type arm is often too small for the contrast every downstream tool assumes.",
        "name": "profile_pdac_driver_landscape",
        "optional_parameters": [
            {
                "default": "Pancreatic Cancer",
                "description": "Cancer context passed to the DepMap lineage matcher",
                "name": "cancer_type",
                "type": "str",
            },
            {
                "default": None,
                "description": "Mutation drivers to profile (list or comma-separated); defaults to the curated "
                "TCGA/ICGC PDAC driver list",
                "name": "mutation_drivers",
                "type": "list[str]",
            },
            {
                "default": None,
                "description": "Amplification drivers to profile; defaults to MYC, GATA6, ERBB2, CCNE1, AKT2",
                "name": "amplification_drivers",
                "type": "list[str]",
            },
            {
                "default": None,
                "description": "Deep-deletion drivers to profile; defaults to CDKN2A, SMAD4, TP53",
                "name": "deletion_drivers",
                "type": "list[str]",
            },
            {
                "default": None,
                "description": "Directory holding the DepMap files (default: ./data/biomni_data/data_lake)",
                "name": "data_lake_path",
                "type": "str",
            },
            {
                "default": None,
                "description": "Optional genotype table overriding the default call source",
                "name": "mutation_csv_path",
                "type": "str",
            },
            {
                "default": None,
                "description": "Write the full landscape table to this CSV path",
                "name": "output_csv_path",
                "type": "str",
            },
        ],
        "required_parameters": [],
    },
    {
        "description": "Map synthetic-lethal candidate genes to compounds that exist, using the Broad Repurposing "
        "Hub, and split the candidate list into targets with a clinical-stage inhibitor, targets with preclinical "
        "tool compounds only, and targets with no compound at all (which require a CRISPR arm rather than a drug "
        "arm). Reports drug names, clinical phase and mechanism of action per gene, with SynLethDB context for "
        "the driver.",
        "name": "map_sl_candidates_to_drugs",
        "optional_parameters": [
            {
                "default": None,
                "description": "Directory holding broad_repurposing_hub_phase_moa_target_info.parquet",
                "name": "data_lake_path",
                "type": "str",
            },
            {
                "default": "Preclinical",
                "description": "Lowest clinical phase to report: 'Launched', 'Phase 3', 'Phase 2', 'Phase 1' or "
                "'Preclinical'",
                "name": "min_clinical_phase",
                "type": "str",
            },
            {
                "default": 4,
                "description": "Number of compounds listed per gene, most clinically advanced first",
                "name": "max_drugs_per_gene",
                "type": "int",
            },
            {
                "default": True,
                "description": "Also report whether each gene is a recorded SynLethDB partner of driver_gene",
                "name": "include_synlethdb_context",
                "type": "bool",
            },
            {
                "default": None,
                "description": "Driver gene for the SynLethDB context lookup, e.g. 'KRAS'",
                "name": "driver_gene",
                "type": "str",
            },
            {
                "default": None,
                "description": "Write the full gene-drug table to this CSV path",
                "name": "output_csv_path",
                "type": "str",
            },
        ],
        "required_parameters": [
            {
                "default": None,
                "description": "Candidate genes (list or comma-separated string), e.g. the CANDIDATE_GENES line "
                "of a discovery run",
                "name": "candidate_genes",
                "type": "list[str]",
            },
        ],
    },
    {
        "description": "Fit measured patient-derived organoid (PDO) dose-response data with a 4-parameter "
        "logistic curve and compare genotype groups. Returns IC50, AUC and Emax per organoid and drug, a Welch "
        "comparison of AUC and IC50 between genotype groups (the comparison that actually tests a synthetic "
        "lethality prediction), QC on failed or censored fits, and an optional dose-response figure.",
        "name": "analyze_organoid_drug_response",
        "optional_parameters": [
            {
                "default": "organoid",
                "description": "Column holding the organoid identifier",
                "name": "organoid_column",
                "type": "str",
            },
            {
                "default": "drug",
                "description": "Column holding the drug name",
                "name": "drug_column",
                "type": "str",
            },
            {
                "default": "concentration_um",
                "description": "Column holding the concentration in uM (must be > 0)",
                "name": "concentration_column",
                "type": "str",
            },
            {
                "default": "viability",
                "description": "Column holding viability as a fraction (0-1) or percent (0-100); the scale is "
                "detected automatically",
                "name": "viability_column",
                "type": "str",
            },
            {
                "default": None,
                "description": "Optional replicate identifier column",
                "name": "replicate_column",
                "type": "str",
            },
            {
                "default": None,
                "description": "CSV mapping organoid to genotype; enables the between-genotype comparison",
                "name": "genotype_csv_path",
                "type": "str",
            },
            {
                "default": "genotype",
                "description": "Genotype column name in the genotype table",
                "name": "genotype_column",
                "type": "str",
            },
            {
                "default": None,
                "description": "Write the per-curve parameter table to this CSV path",
                "name": "output_csv_path",
                "type": "str",
            },
            {
                "default": None,
                "description": "Write a dose-response figure (one panel per drug) to this PNG path",
                "name": "plot_output_path",
                "type": "str",
            },
        ],
        "required_parameters": [
            {
                "default": None,
                "description": "CSV/TSV of the organoid drug screen readout (organoid, drug, concentration, "
                "viability per row)",
                "name": "screen_csv_path",
                "type": "str",
            },
        ],
    },
    {
        "description": "Turn measured organoid CRISPR knockout viability into per-gene effects and a genotype "
        "selectivity test. Each knockout is normalised to the non-targeting control of the same organoid, then "
        "two questions are answered separately: does the knockout reduce viability at all, and is the reduction "
        "larger in driver-mutant organoids than wild-type ones. Emits a per-gene verdict (SELECTIVE, "
        "NON-SELECTIVE core fitness, NO EFFECT, TREND ONLY) and CRISPR_VALIDATION lines for the comparison stage.",
        "name": "analyze_crispr_validation",
        "optional_parameters": [
            {
                "default": "organoid",
                "description": "Column holding the organoid identifier",
                "name": "organoid_column",
                "type": "str",
            },
            {
                "default": "gene",
                "description": "Column holding the targeted gene (or the control label)",
                "name": "gene_column",
                "type": "str",
            },
            {
                "default": "viability",
                "description": "Column holding viability as a fraction or percent",
                "name": "viability_column",
                "type": "str",
            },
            {
                "default": "NTC",
                "description": "Value in the gene column marking the non-targeting control",
                "name": "control_label",
                "type": "str",
            },
            {
                "default": "replicate",
                "description": "Replicate identifier column; replicates are what the per-knockout t-test uses",
                "name": "replicate_column",
                "type": "str",
            },
            {
                "default": None,
                "description": "CSV mapping organoid to genotype; required for the selectivity test",
                "name": "genotype_csv_path",
                "type": "str",
            },
            {
                "default": "genotype",
                "description": "Genotype column name in the genotype table",
                "name": "genotype_column",
                "type": "str",
            },
            {
                "default": None,
                "description": "Which genotype value counts as driver-mutant; by default any label containing "
                "'WT'/'wild' is treated as the control arm",
                "name": "mutant_genotype_label",
                "type": "str",
            },
            {
                "default": None,
                "description": "Write the per-gene selectivity table to this CSV path",
                "name": "output_csv_path",
                "type": "str",
            },
        ],
        "required_parameters": [
            {
                "default": None,
                "description": "CSV/TSV of the organoid CRISPR validation readout (organoid, gene, viability, "
                "replicate per row, including non-targeting control rows)",
                "name": "validation_csv_path",
                "type": "str",
            },
        ],
    },
    {
        "description": "Score Biomni's predicted vulnerabilities against what the organoid experiment measured - "
        "the step that makes the project falsifiable. Joins the predicted ranking to the measured effects and "
        "reports falsified predictions FIRST, then confirmed ones, then experimental hits the prediction missed, "
        "followed by precision@k with a permutation p-value and (only when at least 5 genes are shared) rank "
        "correlation. Ends with an explicit PREDICTIVE / WEAK / NOT PREDICTIVE verdict.",
        "name": "compare_prediction_with_experiment",
        "optional_parameters": [
            {
                "default": "gene",
                "description": "Gene column in the prediction table",
                "name": "prediction_gene_column",
                "type": "str",
            },
            {
                "default": "effect_difference",
                "description": "Score column in the prediction table",
                "name": "prediction_score_column",
                "type": "str",
            },
            {
                "default": "gene",
                "description": "Gene column in the experiment table",
                "name": "experiment_gene_column",
                "type": "str",
            },
            {
                "default": "selectivity_delta",
                "description": "Measured effect column in the experiment table",
                "name": "experiment_effect_column",
                "type": "str",
            },
            {
                "default": True,
                "description": "True when a more negative predicted score means a stronger predicted dependency",
                "name": "prediction_lower_is_stronger",
                "type": "bool",
            },
            {
                "default": True,
                "description": "True when a more negative measured effect means a stronger measured dependency",
                "name": "experiment_lower_is_stronger",
                "type": "bool",
            },
            {
                "default": -0.5,
                "description": "Measured effect at or beyond which a gene counts as an experimental hit; "
                "pre-specify it and do not tune it after seeing the data",
                "name": "experiment_hit_threshold",
                "type": "float",
            },
            {
                "default": 5,
                "description": "Size of the predicted top set used for precision@k",
                "name": "top_k",
                "type": "int",
            },
            {
                "default": 10000,
                "description": "Permutations used to build the precision@k null distribution",
                "name": "n_permutations",
                "type": "int",
            },
            {
                "default": None,
                "description": "Write the joined prediction/experiment table to this CSV path",
                "name": "output_csv_path",
                "type": "str",
            },
        ],
        "required_parameters": [
            {
                "default": None,
                "description": "CSV written by a discovery or ranking tool (gene plus a predicted score column)",
                "name": "prediction_csv_path",
                "type": "str",
            },
            {
                "default": None,
                "description": "CSV written by analyze_crispr_validation, or any table of measured per-gene "
                "effects",
                "name": "experiment_csv_path",
                "type": "str",
            },
        ],
    },
    {
        "description": "Build a self-contained HTML report from the PDAC pipeline outputs: driver testability, "
        "synthetic lethal candidates with their selectivity, druggability classes, live PubMed references per "
        "candidate (real PMIDs with links - genes with no literature are shown as having none), and the organoid "
        "validation / concordance sections once measurements exist. Charts are inline SVG with hover tooltips and "
        "a data table behind each one; no external script or stylesheet is loaded.",
        "name": "generate_pdac_report",
        "optional_parameters": [
            {
                "default": "./pdac_run",
                "description": "Directory holding the pipeline CSVs (driver_landscape.csv, prediction_depmap.csv, "
                "drug_candidates.csv, and optionally pdo_curves.csv, crispr_summary.csv, concordance.csv)",
                "name": "run_dir",
                "type": "str",
            },
            {
                "default": None,
                "description": "Where to write the HTML (default: <run_dir>/pdac_report.html)",
                "name": "output_html_path",
                "type": "str",
            },
            {
                "default": "Pancreatic cancer",
                "description": "Disease context used in the PubMed queries and the report header",
                "name": "disease",
                "type": "str",
            },
            {
                "default": "KRAS",
                "description": "Driver gene used in the PubMed queries and the report header",
                "name": "driver",
                "type": "str",
            },
            {
                "default": True,
                "description": "Query PubMed for per-candidate references (needs internet access)",
                "name": "include_literature",
                "type": "bool",
            },
            {
                "default": 6,
                "description": "PubMed records retrieved per candidate gene",
                "name": "max_papers_per_gene",
                "type": "int",
            },
            {
                "default": 10,
                "description": "Number of candidates charted in detail",
                "name": "top_candidates",
                "type": "int",
            },
            {
                "default": True,
                "description": "True writes a complete HTML document; False writes a body fragment for embedding "
                "in a host that supplies the document skeleton",
                "name": "standalone",
                "type": "bool",
            },
            {
                "default": None,
                "description": "Contact e-mail passed to NCBI Entrez, as recommended by their usage policy",
                "name": "email",
                "type": "str",
            },
            {
                "default": None,
                "description": "NCBI API key; raises the Entrez rate limit",
                "name": "api_key",
                "type": "str",
            },
        ],
        "required_parameters": [],
    },
]
