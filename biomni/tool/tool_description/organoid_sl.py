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
]
