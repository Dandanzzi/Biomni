description = [
    {
        "description": "CALL THIS FIRST, BEFORE committing a synthetic-lethality run to any driver. Report "
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
                "default": True,
                "description": "Include the curated RAS-pathway drug set bundled with the module (KRAS G12C/G12D "
                "inhibitors, SOS1, SHP2, WRN, TEAD and others). The local Repurposing Hub snapshot predates every "
                "KRAS inhibitor, so a KRAS project sees no pharmacology without this",
                "name": "include_curated",
                "type": "bool",
            },
            {
                "default": True,
                "description": "Query ChEMBL live for compounds with a recorded mechanism of action against each "
                "gene; failures degrade to the offline sources and are reported",
                "name": "include_chembl",
                "type": "bool",
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
                "default": None,
                "description": "Scan table (e.g. pancancer_g12d.csv) to run the strict filter cascade over; "
                "defaults to <run_dir>/pancancer_g12d.csv. Adds the protocol funnel, final candidates and the "
                "organoid knockout panel",
                "name": "protocol_csv",
                "type": "str",
            },
            {"default": 0.05, "description": "Cascade gate: BH q-value cutoff", "name": "protocol_q_threshold", "type": "float"},
            {"default": -0.8, "description": "Cascade gate: Cohen's d cutoff", "name": "protocol_min_cohens_d", "type": "float"},
            {"default": 30.0, "description": "Cascade gate: minimum percent of the mutant arm that must be dependent", "name": "protocol_min_pct_mutant", "type": "float"},
            {"default": 50.0, "description": "Cascade gate: maximum percent of all screened lines that may be dependent", "name": "protocol_max_pct_all", "type": "float"},
            {"default": None, "description": "Table from validate_candidates_biologically (default: <run_dir>/biovalidation.csv); adds PDAC-restricted effect, KRAS axis, TCGA survival and structure confidence per candidate", "name": "biovalidation_csv", "type": "str"},
            {"default": 30, "description": "How many near-miss genes to query PubMed for when building the combination-candidate table", "name": "combination_max_genes", "type": "int"},
            {"default": None, "description": "Directory holding the Repurposing Hub table, used for druggable near-miss candidates", "name": "data_lake_path", "type": "str"},
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
    {
        "description": "USE THIS INSTEAD OF a within-cancer-type contrast whenever the driver's wild-type arm "
        "is small (fewer than ~8 lines). Discovers driver-selective dependencies in a pan-cancer cohort, then "
        "checks whether each candidate is actually present in one cancer type. Required when the driver is near-universal in the cancer "
        "of interest and its wild-type arm is too small for a contrast (KRAS in PDAC: 40 mutant vs 4 wild-type, "
        "versus 47 vs 791 pan-cancer for G12D). Supports allele-specific arms, excludes other alleles of the same "
        "driver from both arms, re-tests every candidate with the mutant arm's dominant lineage removed so a "
        "tissue effect is visible, and reports a PRESENT / WEAK / ABSENT IN CONTEXT verdict per gene.",
        "name": "discover_sl_pan_cancer_with_context",
        "optional_parameters": [
            {"default": "KRAS", "description": "Driver gene defining the genotype", "name": "driver_gene", "type": "str"},
            {"default": "G12D", "description": "Restrict the mutant arm to one protein change, e.g. 'G12D'; None uses any mutation", "name": "allele", "type": "str"},
            {"default": "Pancreatic Cancer", "description": "Cancer type the candidates are re-checked in", "name": "context_cancer_type", "type": "str"},
            {"default": None, "description": "Directory holding the DepMap files", "name": "data_lake_path", "type": "str"},
            {"default": None, "description": "Optional genotype table overriding the default call source", "name": "mutation_csv_path", "type": "str"},
            {"default": 20, "description": "Number of candidates carried into the context step and printed", "name": "top_n", "type": "int"},
            {"default": 0.05, "description": "One-sided Welch p-value cutoff", "name": "p_threshold", "type": "float"},
            {"default": 0.25, "description": "Benjamini-Hochberg q-value cutoff", "name": "fdr_threshold", "type": "float"},
            {"default": -0.2, "description": "Required mutant-minus-wild-type gene-effect difference", "name": "min_effect_difference", "type": "float"},
            {"default": -0.3, "description": "The mutant group mean gene effect must be below this", "name": "max_mutant_mean_effect", "type": "float"},
            {"default": None, "description": "Require Cohen's d at or below this value (e.g. -0.8 for a large effect); None disables it", "name": "min_cohens_d", "type": "float"},
            {"default": None, "description": "Require this percentage of the altered arm to be dependent, separating a group-wide shift from one driven by a few lines", "name": "min_pct_mutant_dependent", "type": "float"},
            {"default": None, "description": "Reject genes depleted in more than this percentage of ALL screened lines", "name": "max_pct_all_dependent", "type": "float"},
            {"default": True, "description": "Drop common-essential genes", "name": "exclude_pan_essential", "type": "bool"},
            {"default": None, "description": "Write the full pan-cancer table to this CSV path", "name": "output_csv_path", "type": "str"},
        ],
        "required_parameters": [],
    },
    {
        "description": "Stratify driver-mutant cell lines by a SECOND alteration instead of by the driver itself, "
        "for cancers where the driver is near-universal and has no wild-type control arm. In PDAC this splits "
        "KRAS-mutant lines by TP53, SMAD4 or CDKN2A status and asks what the co-altered genotype depends on. The "
        "result is a dependency of the co-altered genotype, not of the driver, and the report states that "
        "explicitly; both arms carry the driver, so nothing here is evidence about the driver itself.",
        "name": "discover_comutation_stratified_sl",
        "optional_parameters": [
            {"default": "KRAS", "description": "Driver every line in the cohort must carry", "name": "driver_gene", "type": "str"},
            {"default": "TP53", "description": "Second gene whose status splits the cohort, e.g. TP53, SMAD4, CDKN2A", "name": "comutation_gene", "type": "str"},
            {"default": "mutation", "description": "Alteration type for the second gene: mutation, deletion or amplification", "name": "comutation_event", "type": "str"},
            {"default": None, "description": "Restrict the cohort to one driver allele, e.g. 'G12D'", "name": "driver_allele", "type": "str"},
            {"default": "Pancreatic Cancer", "description": "Cancer context", "name": "cancer_type", "type": "str"},
            {"default": None, "description": "Directory holding the DepMap files", "name": "data_lake_path", "type": "str"},
            {"default": None, "description": "Optional genotype table overriding the default call source", "name": "mutation_csv_path", "type": "str"},
            {"default": 20, "description": "Number of candidates printed", "name": "top_n", "type": "int"},
            {"default": 0.05, "description": "One-sided Welch p-value cutoff", "name": "p_threshold", "type": "float"},
            {"default": 0.25, "description": "Benjamini-Hochberg q-value cutoff", "name": "fdr_threshold", "type": "float"},
            {"default": -0.2, "description": "Required between-arm gene-effect difference", "name": "min_effect_difference", "type": "float"},
            {"default": -0.3, "description": "The co-altered group mean gene effect must be below this", "name": "max_mutant_mean_effect", "type": "float"},
            {"default": None, "description": "Require Cohen's d at or below this value (e.g. -0.8 for a large effect); None disables it", "name": "min_cohens_d", "type": "float"},
            {"default": None, "description": "Require this percentage of the altered arm to be dependent, separating a group-wide shift from one driven by a few lines", "name": "min_pct_mutant_dependent", "type": "float"},
            {"default": None, "description": "Reject genes depleted in more than this percentage of ALL screened lines", "name": "max_pct_all_dependent", "type": "float"},
            {"default": True, "description": "Drop common-essential genes", "name": "exclude_pan_essential", "type": "bool"},
            {"default": None, "description": "Write the full table to this CSV path", "name": "output_csv_path", "type": "str"},
        ],
        "required_parameters": [],
    },
    {
        "description": "Separate dependencies shared by two driver alleles from those specific to one, e.g. KRAS "
        "G12D vs G12V. Runs THREE contrasts - allele A vs wild-type, allele B vs wild-type, and allele A vs "
        "allele B directly - and reports a gene as allele-specific only when the direct contrast confirms it, "
        "because a gene can clear the gate for one allele and miss it for the other purely through arm size. "
        "Candidates are checked for presence in a context cancer type and annotated with the MRTX1133 "
        "combination axis they sit on (curated prior knowledge, dated).",
        "name": "compare_allele_specific_dependencies",
        "optional_parameters": [
            {"default": "KRAS", "description": "Driver carrying the alleles", "name": "driver_gene", "type": "str"},
            {"default": "G12D", "description": "First protein change to compare, reported as clinical priority", "name": "allele_a", "type": "str"},
            {"default": "G12V", "description": "Second protein change to compare", "name": "allele_b", "type": "str"},
            {"default": "pan-cancer", "description": "Cohort for the discovery scans: 'pan-cancer' or a cancer type. Within PDAC the wild-type arm is too small for either allele contrast", "name": "discovery_cohort", "type": "str"},
            {"default": "Pancreatic Cancer", "description": "Cancer type candidates are checked for presence in", "name": "context_cancer_type", "type": "str"},
            {"default": None, "description": "Directory holding the DepMap files", "name": "data_lake_path", "type": "str"},
            {"default": None, "description": "Optional genotype table overriding the default call source", "name": "mutation_csv_path", "type": "str"},
            {"default": 15, "description": "Number of rows printed per class", "name": "top_n", "type": "int"},
            {"default": 0.05, "description": "One-sided Welch p-value cutoff", "name": "p_threshold", "type": "float"},
            {"default": 0.25, "description": "Benjamini-Hochberg q-value cutoff", "name": "fdr_threshold", "type": "float"},
            {"default": -0.2, "description": "Required allele-minus-wild-type gene-effect difference", "name": "min_effect_difference", "type": "float"},
            {"default": -0.3, "description": "The allele group mean gene effect must be below this", "name": "max_mutant_mean_effect", "type": "float"},
            {"default": True, "description": "Drop common-essential genes", "name": "exclude_pan_essential", "type": "bool"},
            {"default": None, "description": "Write the joined three-contrast table to this CSV path", "name": "output_csv_path", "type": "str"},
        ],
        "required_parameters": [],
    },
    {
        "description": "Run the biological validation battery on candidate genes: PDAC-restricted effect size "
        "(against both the wild-type arm and all other PDAC lines), KRAS effector axis assignment, TCGA PAAD "
        "survival by median expression split with a log-rank test, direct inhibitor availability OR a druggable "
        "STRING neighbour within one hop, ChEMBL target class, and AlphaFold model confidence (pLDDT - model "
        "confidence and a prerequisite for pocket work, NOT a pocket-quality score).",
        "name": "validate_candidates_biologically",
        "optional_parameters": [
            {"default": "KRAS", "description": "Driver defining the PDAC-restricted contrast", "name": "driver_gene", "type": "str"},
            {"default": "G12D", "description": "Allele defining the PDAC-restricted contrast", "name": "allele", "type": "str"},
            {"default": None, "description": "Directory holding the DepMap files", "name": "data_lake_path", "type": "str"},
            {"default": None, "description": "Optional genotype table overriding the default call source", "name": "mutation_csv_path", "type": "str"},
            {"default": 0.7, "description": "STRING combined-score cutoff for the one-hop neighbourhood", "name": "string_min_score", "type": "float"},
            {"default": 15, "description": "Neighbours examined per gene when looking for a druggable one", "name": "max_neighbours", "type": "int"},
            {"default": True, "description": "Run the TCGA PAAD log-rank test (needs cBioPortal)", "name": "include_survival", "type": "bool"},
            {"default": True, "description": "Fetch AlphaFold model confidence (needs UniProt and AlphaFold DB)", "name": "include_structure", "type": "bool"},
            {"default": True, "description": "Report median GTEx expression in normal pancreas and flag toxicity risk above 50 TPM", "name": "include_gtex", "type": "bool"},
            {"default": True, "description": "Query ClinicalTrials.gov for active Phase 1-3 trials in a pancreatic/KRAS context", "name": "include_trials", "type": "bool"},
            {"default": True, "description": "Replicate the genotype contrast in the independent Sanger Project Score screen", "name": "include_sanger", "type": "bool"},
            {"default": True, "description": "Run fpocket for a binding-pocket druggability score; without structure_paths the AlphaFold model is used, which has no ligands or partners", "name": "include_pockets", "type": "bool"},
            {"default": True, "description": "Test each candidate's drugs and druggable neighbours for allele-selective sensitivity in the PRISM repurposing screen", "name": "include_prism", "type": "bool"},
            {"default": None, "description": "Gene to experimental PDB path, e.g. {'RAB10': '9G0C_RAB10.pdb'}; preferred over AlphaFold when an experimental structure of the right chain exists", "name": "structure_paths", "type": "dict"},
            {"default": None, "description": "Write the per-candidate table to this CSV path", "name": "output_csv_path", "type": "str"},
        ],
        "required_parameters": [
            {"default": None, "description": "Genes to validate (list or comma-separated string)", "name": "candidate_genes", "type": "list[str]"},
        ],
    },
    {
        "description": "Replicate a genotype contrast in the Sanger Project Score CRISPR screen - an independent "
        "library, cell panel and analysis pipeline, so a hit that survives it is not a Broad-specific artefact. "
        "Returns the Sanger effect size per gene with a REPLICATED / NOT REPLICATED verdict against a Cohen's d "
        "threshold. Needs Sanger_ProjectScore_corrected_logFC.parquet in the data lake.",
        "name": "sanger_replication",
        "optional_parameters": [
            {"default": "KRAS", "description": "Driver gene defining the contrast", "name": "driver_gene", "type": "str"},
            {"default": "G12D", "description": "Allele defining the altered arm", "name": "allele", "type": "str"},
            {"default": None, "description": "Directory holding the Sanger matrix", "name": "data_lake_path", "type": "str"},
            {"default": 3, "description": "Minimum lines per arm", "name": "min_group_size", "type": "int"},
            {"default": -0.3, "description": "Cohen's d at or below which the effect counts as replicated", "name": "replication_d_threshold", "type": "float"},
        ],
        "required_parameters": [
            {"default": None, "description": "Genes to test for replication", "name": "genes", "type": "list[str]"},
        ],
    },
    {
        "description": "Audit a filter cascade against known KRAS vulnerabilities (SOS1, EGFR, CDK4, CDK6, PRMT5, "
        "TEAD1, YAP1, RAF1, AURKA, WEE1 by default) to estimate its false-negative rate. Reports the FIRST gate "
        "each known gene fails, its drugs and active trials, and what the gene is established as - a calibration "
        "statement about the gates, not a claim about the genes.",
        "name": "benchmark_known_vulnerabilities",
        "optional_parameters": [
            {"default": None, "description": "Genes to audit; defaults to the curated KRAS benchmark set", "name": "benchmark_genes", "type": "list[str]"},
            {"default": 0.05, "description": "Gate: BH q-value cutoff", "name": "q_threshold", "type": "float"},
            {"default": -0.8, "description": "Gate: Cohen's d cutoff", "name": "min_cohens_d", "type": "float"},
            {"default": 30.0, "description": "Gate: minimum percent of the altered arm that must be dependent", "name": "min_pct_mutant", "type": "float"},
            {"default": 50.0, "description": "Gate: maximum percent of all screened lines that may be dependent", "name": "max_pct_all", "type": "float"},
            {"default": None, "description": "Directory holding the drug tables", "name": "data_lake_path", "type": "str"},
            {"default": False, "description": "Also query ClinicalTrials.gov per gene", "name": "include_trials", "type": "bool"},
            {"default": None, "description": "Write the audit table to this CSV path", "name": "output_csv_path", "type": "str"},
        ],
        "required_parameters": [
            {"default": None, "description": "Scan table from a discovery run (gene, effect_difference, cohens_d, q_value, pct_a_dependent, pct_all_dependent, pan_essential)", "name": "scan_csv_path", "type": "str"},
        ],
    },
    {
        "description": "Test whether a drug preferentially kills driver-allele cell lines in the PRISM "
        "repurposing secondary screen - the pharmacological arm a CRISPR screen cannot supply, since a knockout "
        "says the gene is required while a compound says the molecule works. Lower area under the dose-response "
        "curve means more sensitive, so a negative delta means the allele lines are preferentially killed. Needs "
        "PRISM_secondary_dose_response.csv in the data lake.",
        "name": "prism_allele_sensitivity",
        "optional_parameters": [
            {"default": "KRAS", "description": "Driver gene defining the contrast", "name": "driver_gene", "type": "str"},
            {"default": "G12D", "description": "Allele defining the altered arm", "name": "allele", "type": "str"},
            {"default": None, "description": "Directory holding the PRISM table", "name": "data_lake_path", "type": "str"},
            {"default": None, "description": "Optional genotype table overriding the default call source", "name": "mutation_csv_path", "type": "str"},
            {"default": 3, "description": "Minimum cell lines per arm", "name": "min_group_size", "type": "int"},
        ],
        "required_parameters": [
            {"default": None, "description": "Drug names to test (list or comma-separated string), as they appear in PRISM", "name": "drug_names", "type": "list[str]"},
        ],
    },
]
