description = [
    {
        "description": "Discover candidate synthetic-lethal partner genes of a driver mutation in a given cancer "
        "type by splitting DepMap cell lines into mutant and wild-type groups and testing every gene in the CRISPR "
        "knockout screen with a Welch t-test for genotype-selective dependency. Reports effect sizes, BH-FDR "
        "q-values, selectivity, pan-essentiality QC and full provenance.",
        "name": "discover_synthetic_lethal_candidates",
        "optional_parameters": [
            {
                "default": None,
                "description": "Directory holding DepMap_CRISPRGeneEffect.csv and DepMap_Model.csv "
                "(default: ./data/biomni_data/data_lake)",
                "name": "data_lake_path",
                "type": "str",
            },
            {
                "default": None,
                "description": "Optional mutation table (columns: ModelID, HugoSymbol[, ProteinChange]) that "
                "overrides the default mutation source (local DepMap somatic mutations, else cBioPortal CCLE)",
                "name": "mutation_csv_path",
                "type": "str",
            },
            {
                "default": 20,
                "description": "Number of top candidates to report in detail",
                "name": "top_n",
                "type": "int",
            },
            {
                "default": 0.05,
                "description": "Uncorrected Welch t-test p-value cutoff",
                "name": "p_threshold",
                "type": "float",
            },
            {
                "default": 0.25,
                "description": "Benjamini-Hochberg q-value cutoff; set to 1.0 to disable FDR filtering",
                "name": "fdr_threshold",
                "type": "float",
            },
            {
                "default": -0.2,
                "description": "Required (mutant mean - wild-type mean) gene-effect difference; must be negative",
                "name": "min_effect_difference",
                "type": "float",
            },
            {
                "default": -0.3,
                "description": "The mutant group mean gene effect must be below this value to count as a real "
                "dependency",
                "name": "max_mutant_mean_effect",
                "type": "float",
            },
            {
                "default": True,
                "description": "Drop common-essential genes depleted in more than 80% of all screened cell lines",
                "name": "exclude_pan_essential",
                "type": "bool",
            },
            {
                "default": None,
                "description": "If given, the full ranked candidate table is written to this CSV path",
                "name": "output_csv_path",
                "type": "str",
            },
        ],
        "required_parameters": [
            {
                "default": None,
                "description": "Cancer context to analyse, e.g. 'Pancreatic Cancer', 'Lung', or 'pan-cancer' for "
                "all lineages",
                "name": "cancer_type",
                "type": "str",
            },
            {
                "default": None,
                "description": "HUGO symbol of the mutated driver gene defining the two groups, e.g. 'KRAS'",
                "name": "target_mutation",
                "type": "str",
            },
        ],
    },
    {
        "description": "Validate synthetic-lethal candidate genes against the PubMed literature using the NCBI "
        "Entrez E-utilities API (esearch + efetch with XML abstract parsing). Returns per candidate the query used, "
        "retrieved papers with titles/years/journals/abstracts, and a deterministic 0-100 literature support score "
        "that penalises refuting or non-replication language.",
        "name": "validate_sl_candidates_with_pubmed",
        "optional_parameters": [
            {
                "default": 8,
                "description": "Maximum number of PubMed records to retrieve per candidate gene",
                "name": "max_papers_per_gene",
                "type": "int",
            },
            {
                "default": None,
                "description": "Only consider papers published in or after this year",
                "name": "min_year",
                "type": "int",
            },
            {
                "default": True,
                "description": "Include truncated abstract text for the top papers of each candidate",
                "name": "include_abstracts",
                "type": "bool",
            },
            {
                "default": 700,
                "description": "Number of abstract characters to show per paper",
                "name": "abstract_chars",
                "type": "int",
            },
            {
                "default": None,
                "description": "Contact e-mail passed to Entrez as recommended by NCBI usage policy",
                "name": "email",
                "type": "str",
            },
            {
                "default": None,
                "description": "NCBI API key raising the Entrez rate limit from 3 to 10 requests per second",
                "name": "api_key",
                "type": "str",
            },
        ],
        "required_parameters": [
            {
                "default": None,
                "description": "Disease context used in the query, e.g. 'Pancreatic cancer'",
                "name": "disease",
                "type": "str",
            },
            {
                "default": None,
                "description": "The mutated driver gene, e.g. 'KRAS'",
                "name": "mutated_gene",
                "type": "str",
            },
            {
                "default": None,
                "description": "Candidate synthetic-lethal partner genes; a list of HUGO symbols or a "
                "comma-separated string",
                "name": "candidate_genes",
                "type": "list[str] | str",
            },
        ],
    },
    {
        "description": "Analyse STRING-DB protein-protein interactions between a mutated driver gene and synthetic "
        "lethal candidates. Reports direct edges broken down by evidence channel (experimental, database, "
        "co-expression, text mining), shared first-degree interaction partners as pathway-level proximity, a "
        "proximity classification per candidate, and functional enrichment of the joint network.",
        "name": "analyze_ppi_network_for_sl",
        "optional_parameters": [
            {
                "default": 9606,
                "description": "NCBI taxonomy identifier (9606 = Homo sapiens)",
                "name": "species",
                "type": "int",
            },
            {
                "default": 400,
                "description": "Minimum STRING combined score, 0-1000 (400 = medium confidence)",
                "name": "required_score",
                "type": "int",
            },
            {
                "default": True,
                "description": "Compute shared first-degree neighbours between the driver and each candidate",
                "name": "analyze_shared_partners",
                "type": "bool",
            },
            {
                "default": 250,
                "description": "Maximum interaction partners fetched per gene for the shared-neighbour analysis",
                "name": "partner_limit",
                "type": "int",
            },
            {
                "default": True,
                "description": "Run STRING functional enrichment over the driver plus candidate gene set",
                "name": "run_enrichment",
                "type": "bool",
            },
        ],
        "required_parameters": [
            {
                "default": None,
                "description": "The mutated driver gene, e.g. 'KRAS'",
                "name": "target_gene",
                "type": "str",
            },
            {
                "default": None,
                "description": "Candidate synthetic-lethal partner genes; a list of HUGO symbols or a "
                "comma-separated string",
                "name": "candidate_genes",
                "type": "list[str] | str",
            },
        ],
    },
    {
        "description": "Run the complete synthetic lethality pipeline end to end (DepMap statistical discovery -> "
        "PubMed literature validation -> STRING PPI analysis) and integrate the three evidence streams into an "
        "evidence dossier: per candidate a confidence grade, supporting and contradicting evidence, remaining "
        "uncertainties, a minimal validation experiment and a Go / Hold / No-go recommendation.",
        "name": "generate_sl_evidence_dossier",
        "optional_parameters": [
            {
                "default": 5,
                "description": "Number of top statistical candidates carried into literature and PPI validation",
                "name": "top_n",
                "type": "int",
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
                "default": 8,
                "description": "Maximum PubMed records retrieved per candidate gene",
                "name": "max_papers_per_gene",
                "type": "int",
            },
            {
                "default": 400,
                "description": "Minimum STRING combined score, 0-1000",
                "name": "required_score",
                "type": "int",
            },
            {
                "default": None,
                "description": "Contact e-mail for the NCBI Entrez API",
                "name": "email",
                "type": "str",
            },
            {
                "default": None,
                "description": "NCBI API key for a higher Entrez rate limit",
                "name": "api_key",
                "type": "str",
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
                "description": "Mutated driver gene defining the genotype contrast, e.g. 'KRAS'",
                "name": "target_mutation",
                "type": "str",
            },
        ],
    },
]
