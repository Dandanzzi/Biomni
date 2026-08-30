description = [
    {
        "description": "Rank synthetic lethal partners of a driver gene by combining three independent "
        "inference channels (DAISY design): differential essentiality in driver-deficient lines, "
        "co-expression with the driver, and depletion of double-inactivation among surviving lines. "
        "Returns a RANKED list rather than an FDR-filtered set, because a genome-wide FDR gate "
        "discards true interactions at realistic group sizes. Applies the essentiality floor on DepMap "
        "dependency probabilities rather than gene-effect cutoffs, carries Broad Repurposing Hub "
        "druggability through scoring instead of applying it post hoc, and reports where "
        "already-established interactions rank as a calibration check on the run. Pools multiple driver "
        "genes into one deficiency genotype (e.g. BRCA1+BRCA2) and drops likely-VUS missense calls.",
        "name": "discover_sl_multichannel",
        "optional_parameters": [
            {
                "default": "pan-cancer",
                "description": "Cancer context, or 'pan-cancer'. Pan-cancer is the default because a "
                "single lineage rarely supplies enough driver-mutant lines to power the channels",
                "name": "cancer_type",
                "type": "str",
            },
            {
                "default": None,
                "description": "Interactions already established for this driver (e.g. ['PARP1','POLQ'] "
                "for BRCA); their ranks are reported as a calibration check",
                "name": "known_sl_partners",
                "type": "list[str]",
            },
            {
                "default": True,
                "description": "Count only likely loss-of-function variants as driver-deficient, dropping "
                "simple missense calls that are usually variants of uncertain significance",
                "name": "lof_only",
                "type": "bool",
            },
            {
                "default": 30,
                "description": "Number of ranked candidates to report",
                "name": "top_k",
                "type": "int",
            },
            {
                "default": 3,
                "description": "How many channels a gene must score in to be reported (3 = DAISY's rule)",
                "name": "require_channels",
                "type": "int",
            },
            {
                "default": None,
                "description": "Directory holding the DepMap files",
                "name": "data_lake_path",
                "type": "str",
            },
        ],
        "required_parameters": [
            {
                "default": None,
                "description": "Driver gene(s) whose loss defines the genotype, e.g. ['BRCA1','BRCA2']",
                "name": "driver_genes",
                "type": "list[str]",
            },
        ],
    },
]
