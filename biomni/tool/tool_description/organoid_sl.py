description = [
    {
        "description": "Find synthetic lethal candidates that a pooled cell-line analysis hides by averaging over "
        "transcriptional subtypes. Scores every cell line for subtype from expression (PDAC classical vs "
        "basal-like by default), splits the panel, repeats the mutation-stratified Welch t-test inside each "
        "subtype, and reports genes significant within one subtype but not in the pooled test. These are "
        "interactions a cell-line panel structurally cannot settle, because the relevant subtype is represented "
        "by only a handful of lines - making them the right candidates to take to a subtype-matched "
        "patient-derived organoid.",
        "name": "discover_subtype_masked_sl_candidates",
        "optional_parameters": [
            {
                "default": "pancreatic",
                "description": "Built-in subtype marker set to use: 'pancreatic' (classical vs basal-like) or "
                "'colorectal'",
                "name": "subtype_marker_set",
                "type": "str",
            },
            {
                "default": None,
                "description": "Directory holding the DepMap files (default: ./data/biomni_data/data_lake)",
                "name": "data_lake_path",
                "type": "str",
            },
            {
                "default": None,
                "description": "Optional mutation table overriding the default mutation source",
                "name": "mutation_csv_path",
                "type": "str",
            },
            {
                "default": 0.05,
                "description": "Uncorrected Welch p-value cutoff applied inside each subtype",
                "name": "p_threshold",
                "type": "float",
            },
            {
                "default": -0.2,
                "description": "Required mutant-minus-wild-type gene-effect difference; must be negative",
                "name": "min_effect_difference",
                "type": "float",
            },
            {
                "default": -0.3,
                "description": "The mutant group mean gene effect must be below this value",
                "name": "max_mutant_mean_effect",
                "type": "float",
            },
            {
                "default": 15,
                "description": "Number of subtype-masked candidates to report",
                "name": "top_n",
                "type": "int",
            },
        ],
        "required_parameters": [
            {
                "default": None,
                "description": "Cancer context, e.g. 'Pancreatic Cancer'",
                "name": "cancer_type",
                "type": "str",
            },
            {
                "default": None,
                "description": "Driver gene defining the genotype contrast, e.g. 'KRAS'",
                "name": "target_mutation",
                "type": "str",
            },
        ],
    },
    {
        "description": "Predict whether a synthetic lethal hit found in 2D cell lines will still be detectable in "
        "a patient-derived organoid, before committing to the experiment. Runs four checks per candidate: niche "
        "coupling against the factors standard organoid medium supplies (WNT3A/RSPO1, EGF, FGF10, Noggin, A83-01, "
        "Y-27632) which can rescue a knockout and cause a false negative; anchorage sensitivity measured from "
        "adherent versus suspension cell lines; therapeutic window from GTEx normal-tissue expression; and "
        "orthogonal human genetic interactions from BioGRID. Returns an ORGANOID-ENHANCED / ORGANOID-MASKED / "
        "ORGANOID-WEAKENED / TRANSFERABLE / UNCERTAIN verdict plus the medium modification needed to keep the "
        "experiment interpretable.",
        "name": "assess_organoid_transferability",
        "optional_parameters": [
            {
                "default": "KRAS",
                "description": "The driver gene the candidates were called against",
                "name": "target_mutation",
                "type": "str",
            },
            {
                "default": "Pancreatic Cancer",
                "description": "Cancer context, used for reporting",
                "name": "cancer_type",
                "type": "str",
            },
            {
                "default": "Pancreas",
                "description": "GTEx tissue used as the normal counterpart for the therapeutic-window check",
                "name": "normal_tissue",
                "type": "str",
            },
            {
                "default": None,
                "description": "Directory holding the DepMap and GTEx files",
                "name": "data_lake_path",
                "type": "str",
            },
        ],
        "required_parameters": [
            {
                "default": None,
                "description": "Candidate genes to assess; a list of HUGO symbols or a comma-separated string",
                "name": "candidate_genes",
                "type": "list[str] | str",
            },
        ],
    },
    {
        "description": "Find dependencies that a serum-grown cell-line CRISPR screen cannot see but a "
        "serum-free organoid screen can. The Nature 2026 tumour-derived organoid biobank found 97 "
        "organoid-specific core fitness genes enriched for steroid/cholesterol/isoprenoid biosynthesis; the "
        "mechanism is that organoid medium is serum-free, so sterols and fatty acids must be made de novo, "
        "while a screen in 10% FBS supplies them. DepMap annotates which lines are cultured serum-free, so "
        "the contrast is measured rather than assumed: each serum-free line is z-scored against "
        "serum-cultured lines of its own lineage and growth pattern. Reports a permutation test per lipid "
        "programme (with lipoprotein uptake as an inverse control), then ranks individual genes and flags the "
        "ones that are NOT dependencies in the 2D cancer panel - the candidates a cell-line screen of this "
        "tumour type would have discarded.",
        "name": "discover_serum_masked_dependencies",
        "optional_parameters": [
            {
                "default": "Pancreatic Cancer",
                "description": "Cancer context whose 2D dependency is the 'would we have found it?' reference",
                "name": "cancer_type",
                "type": "str",
            },
            {
                "default": "KRAS",
                "description": "Driver gene; the mutant-vs-wild-type difference is reported alongside each "
                "candidate. Pass None to skip genotype stratification",
                "name": "target_mutation",
                "type": "str",
            },
            {
                "default": None,
                "description": "Directory holding the DepMap files (default: ./data/biomni_data/data_lake)",
                "name": "data_lake_path",
                "type": "str",
            },
            {
                "default": None,
                "description": "Optional mutation table overriding the default mutation source",
                "name": "mutation_csv_path",
                "type": "str",
            },
            {
                "default": 5000,
                "description": "Permutations used for the programme-level enrichment test",
                "name": "n_permutations",
                "type": "int",
            },
            {
                "default": 4,
                "description": "A gene must be measured in at least this many serum-free lines to be ranked",
                "name": "min_models_scored",
                "type": "int",
            },
            {
                "default": -1.0,
                "description": "Lineage-matched z-score at or below which a gene counts as serum-masked",
                "name": "z_threshold",
                "type": "float",
            },
            {
                "default": -0.5,
                "description": "The gene must also be a real dependency in the serum-free lines: mean gene "
                "effect there at or below this value",
                "name": "max_serum_free_effect",
                "type": "float",
            },
            {
                "default": 1.0,
                "description": "Minimum mean expression in the cancer cohort, log2(TPM+1), to exclude "
                "copy-number and multi-mapping artefacts",
                "name": "min_cohort_expression",
                "type": "float",
            },
            {
                "default": 20,
                "description": "Number of serum-masked candidates to report",
                "name": "top_n",
                "type": "int",
            },
            {
                "default": 80.0,
                "description": "Drop genes essential in more than this percentage of all screened lines",
                "name": "max_pct_all_dependent",
                "type": "float",
            },
            {
                "default": None,
                "description": "Write the full ranked candidate table to this CSV path as well",
                "name": "output_csv_path",
                "type": "str",
            },
        ],
        "required_parameters": [],
    },
    {
        "description": "Find synthetic lethal candidates that the label '<gene>-mutant' hides by pooling "
        "alleles. The Nature 2026 organoid biobank showed KRAS G12 organoids depend on KRAS, EGFR and PTPN11 "
        "while KRAS Q61H organoids ignore EGFR inhibition and EGF withdrawal - a split a pooled "
        "mutant-vs-wild-type test averages away. Splits the mutant lines by protein change and runs two "
        "contrasts per allele: allele vs wild type, and allele vs the other mutant alleles (the direct test of "
        "what pooling hides, which needs no wild-type lines). Also reports a per-allele read-out of RAS "
        "pathway genes, and lists the alleles present in patients but too rare in the cell-line panel to test - "
        "the specific gap a patient-derived organoid biobank can be built to fill.",
        "name": "discover_allele_resolved_sl_candidates",
        "optional_parameters": [
            {
                "default": "Pancreatic Cancer",
                "description": "Cancer context, e.g. 'Pancreatic Cancer'",
                "name": "cancer_type",
                "type": "str",
            },
            {
                "default": "KRAS",
                "description": "Driver gene whose alleles define the groups",
                "name": "target_mutation",
                "type": "str",
            },
            {
                "default": None,
                "description": "Directory holding the DepMap files (default: ./data/biomni_data/data_lake)",
                "name": "data_lake_path",
                "type": "str",
            },
            {
                "default": None,
                "description": "Optional mutation table overriding the default source; needs a ProteinChange "
                "column for allele resolution",
                "name": "mutation_csv_path",
                "type": "str",
            },
            {
                "default": 4,
                "description": "Minimum cell lines carrying an allele before it is tested",
                "name": "min_allele_lines",
                "type": "int",
            },
            {
                "default": 0.05,
                "description": "Uncorrected Welch p-value cutoff",
                "name": "p_threshold",
                "type": "float",
            },
            {
                "default": -0.2,
                "description": "Required gene-effect difference against the comparison group; must be negative",
                "name": "min_effect_difference",
                "type": "float",
            },
            {
                "default": -0.3,
                "description": "The allele group mean gene effect must be below this value",
                "name": "max_mutant_mean_effect",
                "type": "float",
            },
            {
                "default": 80.0,
                "description": "Drop genes essential in more than this percentage of all screened lines",
                "name": "max_pct_all_dependent",
                "type": "float",
            },
            {
                "default": 15,
                "description": "Number of allele-restricted candidates to report",
                "name": "top_n",
                "type": "int",
            },
            {
                "default": None,
                "description": "Write the full allele-resolved result table to this CSV path as well",
                "name": "output_csv_path",
                "type": "str",
            },
        ],
        "required_parameters": [],
    },
    {
        "description": "Write a patient-derived organoid (PDO) protocol that can confirm or kill a synthetic "
        "lethal candidate: model panel with subtype coverage and a matched normal-organoid arm, medium "
        "formulation adjusted for the candidate's niche coupling, CRISPR delivery into organoids, 3D-appropriate "
        "readouts, powering with the patient line as the unit of replication, and pre-specified falsification "
        "criteria. Also states explicitly what the experiment cannot settle (stroma, immune contribution, "
        "pharmacology).",
        "name": "design_organoid_sl_experiment",
        "optional_parameters": [
            {
                "default": "KRAS",
                "description": "The driver gene defining the genotype contrast",
                "name": "target_mutation",
                "type": "str",
            },
            {
                "default": "Pancreatic Cancer",
                "description": "Cancer context",
                "name": "cancer_type",
                "type": "str",
            },
            {
                "default": "both",
                "description": "Which transcriptional subtype the organoid panel should cover: 'classical', "
                "'basal' or 'both'",
                "name": "subtype",
                "type": "str",
            },
            {
                "default": None,
                "description": "Directory holding the DepMap files, used to list organoid models present in DepMap",
                "name": "data_lake_path",
                "type": "str",
            },
        ],
        "required_parameters": [
            {
                "default": None,
                "description": "The candidate synthetic lethal partner to test, e.g. 'TEAD1'",
                "name": "candidate_gene",
                "type": "str",
            },
        ],
    },
    {
        "description": "Rank organoid synthetic-lethal candidates against one another and return a single "
        "recommendation of which candidate an organoid experiment should be spent on. The per-axis tools "
        "each emit an independent verdict per gene and never compare candidates, so the decision had to be "
        "made by hand outside the tooling; this makes it a deterministic rule. Two distinctions drive the "
        "ranking and neither is made elsewhere: (1) genotype selectivity - a dependency equally strong in "
        "driver-wild-type cells is organoid-specific core fitness, not synthetic lethality, and the serum "
        "axis selects on the medium effect alone so such genes reach the candidate list looking like "
        "discoveries; (2) why the organoid helps - ORGANOID-ENHANCED is awarded both for anchorage in matrix "
        "(the 3D format measures the same question better) and for the absence of serum (the medium changes "
        "which metabolic genes are limiting, a window risk rather than a target). Also applies a graded "
        "paralog check, a normal-tissue therapeutic-window check and a bounded literature modulator, and "
        "prints every threshold it used.",
        "name": "rank_organoid_sl_candidates",
        "optional_parameters": [
            {
                "default": "KRAS",
                "description": "The driver gene defining the genotype contrast",
                "name": "target_mutation",
                "type": "str",
            },
            {
                "default": "Pancreatic Cancer",
                "description": "Cancer context used to select the cell-line cohort",
                "name": "cancer_type",
                "type": "str",
            },
            {
                "default": None,
                "description": "Mapping of gene to the discovery axes it came from, e.g. {'SCAP': ['allele', "
                "'lipid']}; convergence across independent axes earns a bounded bonus",
                "name": "axis_origins",
                "type": "dict",
            },
            {
                "default": None,
                "description": "Mapping of gene to the 0-100 score from validate_sl_candidates_with_pubmed; "
                "used only as a modulator, and absence is treated as neutral rather than as evidence against",
                "name": "literature_scores",
                "type": "dict",
            },
            {
                "default": None,
                "description": "Mapping of gene to the genome-wide FDR q-value from "
                "discover_synthetic_lethal_candidates; omitted rather than recomputed, because an FDR over a "
                "handful of candidates would not mean what a genome-wide one does",
                "name": "q_values",
                "type": "dict",
            },
            {
                "default": "Pancreas",
                "description": "GTEx tissue used as the normal counterpart for the therapeutic-window check",
                "name": "normal_tissue",
                "type": "str",
            },
            {
                "default": None,
                "description": "Directory holding the DepMap and GTEx files",
                "name": "data_lake_path",
                "type": "str",
            },
            {
                "default": None,
                "description": "Optional mutation table overriding the default mutation source",
                "name": "mutation_csv_path",
                "type": "str",
            },
            {
                "default": 25,
                "description": "Number of ranked rows to print",
                "name": "top_n",
                "type": "int",
            },
        ],
        "required_parameters": [
            {
                "default": None,
                "description": "Candidates to rank; a list, or a comma/whitespace separated string",
                "name": "candidate_genes",
                "type": "list[str]",
            },
        ],
    },
]
