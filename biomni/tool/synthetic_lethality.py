"""Synthetic lethality (SL) discovery and evidence-validation tools for Biomni.

This module implements a three-stage tool chain for context-specific synthetic
lethality research (e.g. "KRAS-mutant pancreatic cancer"):

1. ``discover_synthetic_lethal_candidates`` - statistical discovery of candidate
   SL partners from DepMap CRISPR gene-effect data stratified by mutation status.
2. ``validate_sl_candidates_with_pubmed`` - literature validation of candidates via
   the NCBI Entrez E-utilities API (esearch + efetch abstract parsing).
3. ``analyze_ppi_network_for_sl`` - protein-protein interaction validation via the
   STRING-DB REST API.
4. ``check_dependency_confounders`` - tests whether a candidate's dependency is really
   explained by the driver mutation or by a confounder (paralog loss, lineage, expression
   biomarker), before any wet-lab commitment.
5. ``generate_sl_evidence_dossier`` - orchestrates 1-4 into a single evidence dossier
   with Go / Hold / No-go recommendations.

Design notes
------------
* All numeric judgements are made by explicit code (statistics, thresholds, QC),
  never by free-text summarisation, so the LLM layer only plans and narrates.
* Every tool reports provenance (data file, release, sample counts, API endpoint)
  and QC warnings (pan-essentiality, low sample size, lineage confounding) so that
  weak candidates can be actively falsified rather than silently promoted.
"""

import json
import os
import re
import time
import urllib.parse
import xml.etree.ElementTree as ET
from datetime import datetime, timezone

import requests

# Aliased so the module also imports on Python 3.10; datetime.UTC only exists from 3.11.
_UTC = timezone.utc  # noqa: UP017

DEFAULT_DATA_LAKE = "./data/biomni_data/data_lake"
CBIOPORTAL_API = "https://www.cbioportal.org/api"
CCLE_STUDY_ID = "ccle_broad_2019"
STRING_API = "https://string-db.org/api"
ENTREZ_API = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"

# Cell lines whose gene effect is below this value are treated as "depleted".
DEPLETION_THRESHOLD = -0.5
# A gene depleted in more than this fraction of ALL screened lines is pan-essential.
PAN_ESSENTIAL_FRACTION = 0.80

# Module level caches so that the ~400 MB DepMap matrices are read only once per session.
_DEPMAP_CACHE: dict = {}
_MUTATION_CACHE: dict = {}
_EXPRESSION_CACHE: dict = {}


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------
def _normalize_cell_line_name(name: str) -> str:
    """Normalise a cell line name for cross-resource matching (MIA PaCa-2 -> MIAPACA2)."""
    return re.sub(r"[^A-Z0-9]", "", str(name).upper())


def _http_get(url: str, params: dict | None = None, timeout: int = 60, retries: int = 3, as_json: bool = True):
    """GET with exponential backoff. Returns parsed JSON, raw text, or None on failure."""
    last_error = None
    for attempt in range(retries):
        try:
            response = requests.get(url, params=params, timeout=timeout)
            response.raise_for_status()
            return response.json() if as_json else response.text
        except Exception as e:  # network layer: any failure is retried
            last_error = e
            time.sleep(1.5 * (attempt + 1))
    print(f"[synthetic_lethality] request failed after {retries} attempts ({url}): {last_error}")
    return None


def _resolve_data_lake(data_lake_path: str | None) -> str:
    """Return a usable data lake directory containing the DepMap files."""
    candidates = [
        data_lake_path,
        DEFAULT_DATA_LAKE,
        os.path.join(os.path.expanduser("~"), "biomni_data", "data_lake"),
        "./biomni_data/data_lake",
    ]
    for candidate in candidates:
        if candidate and os.path.exists(os.path.join(candidate, "DepMap_CRISPRGeneEffect.csv")):
            return candidate
    raise FileNotFoundError(
        "DepMap_CRISPRGeneEffect.csv / DepMap_Model.csv were not found. Provide `data_lake_path` "
        f"explicitly (searched: {[c for c in candidates if c]})."
    )


def _file_provenance(path: str) -> str:
    stat = os.stat(path)
    modified = datetime.fromtimestamp(stat.st_mtime, tz=_UTC).strftime("%Y-%m-%d")
    return f"{os.path.basename(path)} ({stat.st_size / 1e6:.0f} MB, local snapshot dated {modified})"


def _load_depmap(data_lake_path: str | None = None) -> dict:
    """Load and cache the DepMap CRISPR gene-effect matrix and the model annotation table.

    Returns
    -------
    dict with keys ``gene_effect`` (DataFrame: ModelID x gene), ``model`` (DataFrame),
    ``gene_columns`` (mapping HUGO symbol -> matrix column) and ``provenance`` (list of str).
    """
    import pandas as pd

    resolved = _resolve_data_lake(data_lake_path)
    if resolved in _DEPMAP_CACHE:
        return _DEPMAP_CACHE[resolved]

    effect_path = os.path.join(resolved, "DepMap_CRISPRGeneEffect.csv")
    model_path = os.path.join(resolved, "DepMap_Model.csv")

    gene_effect = pd.read_csv(effect_path, index_col=0)
    gene_effect.index.name = "ModelID"
    model = pd.read_csv(model_path)

    # DepMap columns look like "KRAS (3845)"; map the HUGO symbol to the full column name.
    gene_columns = {}
    for column in gene_effect.columns:
        gene_columns[column.split(" (")[0].strip().upper()] = column

    bundle = {
        "gene_effect": gene_effect,
        "model": model,
        "gene_columns": gene_columns,
        "data_lake_path": resolved,
        "provenance": [
            f"CRISPR gene effect: {_file_provenance(effect_path)}; "
            f"{gene_effect.shape[0]} models x {gene_effect.shape[1]} genes (Chronos-corrected)",
            f"Model annotation: {_file_provenance(model_path)}; {model.shape[0]} models",
        ],
    }
    _DEPMAP_CACHE[resolved] = bundle
    return bundle


def _load_expression(data_lake_path: str | None = None) -> dict:
    """Load and cache the DepMap protein-coding expression matrix (log2 TPM+1)."""
    import pandas as pd

    resolved = _resolve_data_lake(data_lake_path)
    if resolved in _EXPRESSION_CACHE:
        return _EXPRESSION_CACHE[resolved]

    path = os.path.join(resolved, "DepMap_OmicsExpressionProteinCodingGenesTPMLogp1.csv")
    if not os.path.exists(path):
        raise FileNotFoundError(f"Expression matrix not found at {path}; confounder analysis needs it.")

    expression = pd.read_csv(path, index_col=0)
    expression.index.name = "ModelID"
    bundle = {
        "expression": expression,
        "gene_columns": {c.split(" (")[0].strip().upper(): c for c in expression.columns},
        "provenance": (
            f"Expression: {_file_provenance(path)}; {expression.shape[0]} models x {expression.shape[1]} genes "
            "(log2 TPM+1)"
        ),
    }
    _EXPRESSION_CACHE[resolved] = bundle
    return bundle


# SynLethDB 2.0 (synlethdb_human_sl.parquet: gene_a, gene_b, score, source).
# The `score` column of the local snapshot is constant 1.0 for all 37,943 pairs, so it carries
# no confidence information; the evidence type in `source` is the only usable ranking signal.
_SYNLETHDB_EVIDENCE_TIER = {
    "LOW THROUGHPUT": 3,
    "CRISPR/CRISPRI": 3,
    "DRUG SCREEN": 3,
    "SYNLETHALITY": 3,
    "HIGH THROUGHPUT": 2,
    "GENOMERNAI": 2,
    "RNAI SCREEN": 2,
    "TEXT MINING": 1,
    "DECIPHER": 1,
    "DAISY": 0,
    "COMPUTATIONAL PREDICTION": 0,
}
_SYNLETHDB_TIER_NAME = {
    3: "experimental (low-throughput / CRISPR / drug screen)",
    2: "high-throughput or RNAi screen",
    1: "text mining / curation",
    0: "computational prediction",
}
_SYNLETHDB_CACHE: dict = {}


def _load_synlethdb(data_lake_path: str | None = None) -> dict | None:
    """Load and cache the human SynLethDB SL pair table. Returns None when the file is absent."""
    import pandas as pd

    try:
        resolved = _resolve_data_lake(data_lake_path)
    except FileNotFoundError:
        resolved = data_lake_path or DEFAULT_DATA_LAKE
    if resolved in _SYNLETHDB_CACHE:
        return _SYNLETHDB_CACHE[resolved]

    path = os.path.join(resolved, "synlethdb_human_sl.parquet")
    if not os.path.exists(path):
        _SYNLETHDB_CACHE[resolved] = None
        return None

    table = pd.read_parquet(path)
    missing = {"gene_a", "gene_b"} - set(table.columns)
    if missing:
        raise ValueError(f"{path} is missing required column(s) {sorted(missing)}")
    table["gene_a"] = table["gene_a"].astype(str).str.upper().str.strip()
    table["gene_b"] = table["gene_b"].astype(str).str.upper().str.strip()
    if "source" not in table.columns:
        table["source"] = ""
    constant_score = "score" in table.columns and table["score"].nunique(dropna=True) <= 1

    bundle = {
        "table": table,
        "constant_score": constant_score,
        "provenance": (
            f"SynLethDB human SL pairs: {_file_provenance(path)}; {len(table)} pairs, "
            f"{len(set(table['gene_a']) | set(table['gene_b']))} genes"
            + (
                " (WARNING: the `score` column is constant in this snapshot, so pairs are weighted by "
                "evidence type instead)"
                if constant_score
                else ""
            )
        ),
    }
    _SYNLETHDB_CACHE[resolved] = bundle
    return bundle


def _synlethdb_partners(gene: str, data_lake_path: str | None = None) -> dict:
    """Return ``{partner: {"sources": [...], "tier": int}}`` for known SL partners of ``gene``.

    An empty dict means either "no partners recorded" or "SynLethDB not available"; callers should
    check :func:`_load_synlethdb` separately when they need to tell those apart.
    """
    bundle = _load_synlethdb(data_lake_path)
    if bundle is None:
        return {}

    gene = str(gene).strip().upper()
    table = bundle["table"]
    hits = table[(table["gene_a"] == gene) | (table["gene_b"] == gene)]

    partners: dict[str, dict] = {}
    for _, row in hits.iterrows():
        partner = row["gene_b"] if row["gene_a"] == gene else row["gene_a"]
        if partner == gene or partner in {"", "NAN", "NONE"}:
            continue
        # SynLethDB concatenates several evidence types with ';' or '|'.
        sources = [s.strip() for s in re.split(r"[;|]", str(row["source"])) if s.strip()]
        tier = max((_SYNLETHDB_EVIDENCE_TIER.get(s.upper(), 0) for s in sources), default=0)
        record = partners.setdefault(partner, {"sources": set(), "tier": 0})
        record["sources"].update(sources)
        record["tier"] = max(record["tier"], tier)
    for record in partners.values():
        record["sources"] = sorted(record["sources"])
    return partners


_OVARIAN_HISTOLOGY_ALIASES = {
    "hgsoc": "High-Grade Serous",
    "hgsc": "High-Grade Serous",
    "high grade serous": "High-Grade Serous",
    "high-grade serous": "High-Grade Serous",
    "serous": "Serous",
    "clear cell": "Clear Cell",
    "ccoc": "Clear Cell",
    "endometrioid": "Endometrioid",
    "mucinous": "Mucinous",
    "germ cell": "Germ Cell",
}


def _select_ovarian_models(model_df, histology: str | None = None):
    """Subset the DepMap model table to ovarian cancer lines, optionally to one histology.

    Ovarian lines sit under OncotreeLineage "Ovary/Fallopian Tube"; the histology (high-grade
    serous, clear cell, endometrioid, mucinous) lives in OncotreeSubtype. Non-cancerous and
    immortalized normal ovarian lines are always dropped - they have no driver genotype to
    stratify on and would inflate the wild-type group.
    """
    mask = model_df["OncotreeLineage"].astype(str).str.contains("ovar", case=False, na=False)
    for column in ("OncotreePrimaryDisease", "OncotreeSubtype"):
        if column in model_df.columns:
            mask = mask | model_df[column].astype(str).str.contains("ovarian", case=False, na=False)

    subset = model_df[mask].copy()
    non_cancerous = subset["OncotreePrimaryDisease"].astype(str).str.contains(
        "non-cancerous", case=False, na=False
    ) | subset["OncotreeSubtype"].astype(str).str.contains("immortalized", case=False, na=False)
    n_dropped = int(non_cancerous.sum())
    subset = subset[~non_cancerous]
    note = f"OncotreeLineage 'Ovary/Fallopian Tube' -> {len(subset)} tumour lines ({n_dropped} non-cancerous dropped)"

    if histology:
        key = _OVARIAN_HISTOLOGY_ALIASES.get(str(histology).strip().lower(), str(histology).strip())
        histology_mask = subset["OncotreeSubtype"].astype(str).str.contains(key, case=False, na=False)
        subset = subset[histology_mask]
        note += f"; histology filter '{histology}' -> '{key}' -> {len(subset)} lines"
    return subset, note


def _select_cancer_models(model_df, cancer_type: str):
    """Subset the DepMap model table to a cancer context.

    Matches ``cancer_type`` as a case-insensitive substring against OncotreeLineage,
    OncotreePrimaryDisease, OncotreeSubtype and OncotreeCode. ``"pan-cancer"`` / ``"all"``
    returns every model.
    """
    if str(cancer_type).strip().lower() in {"pan-cancer", "pancancer", "all", "any"}:
        return model_df.copy(), "pan-cancer (no lineage filter)"

    query = str(cancer_type).strip().lower()
    # "pancreatic cancer" -> "pancrea" so that lineage/disease/subtype spellings all match.
    stem = re.sub(r"\s*(cancer|carcinoma|tumou?r|adenocarcinoma|neoplasm)s?\s*$", "", query).strip()
    stem = stem[:7] if len(stem) > 7 else stem

    columns = ["OncotreeLineage", "OncotreePrimaryDisease", "OncotreeSubtype", "OncotreeCode"]
    mask = None
    for column in columns:
        if column not in model_df.columns:
            continue
        hit = model_df[column].astype(str).str.contains(stem, case=False, na=False, regex=False)
        mask = hit if mask is None else (mask | hit)

    subset = model_df[mask].copy() if mask is not None else model_df.iloc[0:0].copy()
    lineages = sorted(subset["OncotreeLineage"].dropna().unique().tolist()) if len(subset) else []
    return subset, f"matched '{cancer_type}' (stem '{stem}') -> lineages {lineages}"


def _mutation_status_from_local_file(gene: str, data_lake_path: str) -> dict | None:
    """Read mutation status from a local DepMap OmicsSomaticMutations.csv if the user has one."""
    import pandas as pd

    for filename in ("DepMap_OmicsSomaticMutations.csv", "OmicsSomaticMutations.csv"):
        path = os.path.join(data_lake_path, filename)
        if not os.path.exists(path):
            continue
        mutations = pd.read_csv(path, usecols=lambda c: c in {"ModelID", "HugoSymbol", "ProteinChange", "VariantInfo"})
        hits = mutations[mutations["HugoSymbol"].astype(str).str.upper() == gene.upper()]
        variants = {}
        for _, row in hits.iterrows():
            variants.setdefault(row["ModelID"], []).append(str(row.get("ProteinChange", "NA")))
        return {
            "mutant_variants": {k: ",".join(sorted(set(v))) for k, v in variants.items()},
            "profiled_models": set(mutations["ModelID"].unique()),
            "key": "ModelID",
            "source": f"local {filename} ({_file_provenance(path)})",
        }
    return None


def _mutation_status_from_cbioportal(gene: str) -> dict | None:
    """Fetch cell-line mutation calls for ``gene`` from the cBioPortal CCLE study.

    Returns mutant variants and the set of profiled samples keyed by normalised cell line
    name, so that "profiled but not mutated" (= wild type) can be distinguished from
    "never profiled" (= unknown).
    """
    cache_key = gene.upper()
    if cache_key in _MUTATION_CACHE:
        return _MUTATION_CACHE[cache_key]

    gene_info = _http_get(f"{CBIOPORTAL_API}/genes/{urllib.parse.quote(gene.upper())}")
    if not gene_info or "entrezGeneId" not in gene_info:
        return None

    mutations = _http_get(
        f"{CBIOPORTAL_API}/molecular-profiles/{CCLE_STUDY_ID}_mutations/mutations",
        params={
            "sampleListId": f"{CCLE_STUDY_ID}_all",
            "entrezGeneId": gene_info["entrezGeneId"],
            "projection": "SUMMARY",
        },
        timeout=90,
    )
    samples = _http_get(f"{CBIOPORTAL_API}/studies/{CCLE_STUDY_ID}/samples", timeout=90)
    if mutations is None or samples is None:
        return None

    variants: dict[str, list[str]] = {}
    for record in mutations:
        key = _normalize_cell_line_name(str(record.get("sampleId", "")).split("_")[0])
        change = record.get("proteinChange") or record.get("mutationType") or "mutated"
        variants.setdefault(key, []).append(str(change))

    result = {
        "mutant_variants": {k: ",".join(sorted(set(v))) for k, v in variants.items()},
        "profiled_models": {_normalize_cell_line_name(str(s["sampleId"]).split("_")[0]) for s in samples},
        "key": "normalized_name",
        "source": (
            f"cBioPortal {CCLE_STUDY_ID} mutation calls "
            f"({len(mutations)} {gene.upper()} variants across {len(samples)} profiled cell lines)"
        ),
    }
    _MUTATION_CACHE[cache_key] = result
    return result


def _annotate_mutation_status(models, gene: str, data_lake_path: str, mutation_csv_path: str | None = None):
    """Attach MUT / WT / UNKNOWN status for ``gene`` to a DepMap model subset.

    Precedence: user-supplied CSV > local DepMap somatic mutation file > cBioPortal CCLE API.
    """
    import pandas as pd

    models = models.copy()
    models["_norm_name"] = models["StrippedCellLineName"].map(_normalize_cell_line_name)

    annotation = None
    if mutation_csv_path:
        table = pd.read_csv(mutation_csv_path)
        required = {"ModelID", "HugoSymbol"}
        if not required.issubset(table.columns):
            raise ValueError(f"{mutation_csv_path} must contain columns {sorted(required)}")
        hits = table[table["HugoSymbol"].astype(str).str.upper() == gene.upper()]
        annotation = {
            "mutant_variants": {r["ModelID"]: str(r.get("ProteinChange", "mutated")) for _, r in hits.iterrows()},
            "profiled_models": set(table["ModelID"].unique()),
            "key": "ModelID",
            "source": f"user-supplied mutation table {mutation_csv_path}",
        }
    if annotation is None:
        annotation = _mutation_status_from_local_file(gene, data_lake_path)
    if annotation is None:
        annotation = _mutation_status_from_cbioportal(gene)
    if annotation is None:
        raise RuntimeError(
            f"Could not determine {gene} mutation status: no local mutation table and the cBioPortal "
            "API is unreachable. Supply `mutation_csv_path` with columns ModelID,HugoSymbol[,ProteinChange]."
        )

    key_column = "ModelID" if annotation["key"] == "ModelID" else "_norm_name"

    def classify(key):
        if key in annotation["mutant_variants"]:
            return "MUT"
        return "WT" if key in annotation["profiled_models"] else "UNKNOWN"

    models["MutationStatus"] = models[key_column].map(classify)
    models["Variant"] = models[key_column].map(lambda k: annotation["mutant_variants"].get(k, ""))
    return models, annotation["source"]


# cBioPortal discrete copy-number calls (GISTIC-style): 2 = amplification, -2 = deep deletion.
CN_AMPLIFICATION = 2
CN_DEEP_DELETION = -2
# DepMap OmicsCNGene stores log2(relative copy number + 1), so a neutral diploid gene is 1.0.
# Amplified = relative CN >= 2 (log2(3) = 1.585); deep deletion = relative CN <= 0.25 (log2(1.25) = 0.322).
DEPMAP_CN_AMPLIFIED = 1.585
DEPMAP_CN_DELETED = 0.322

_CN_CACHE: dict = {}


def _http_post(url: str, payload: dict, params: dict | None = None, timeout: int = 90, retries: int = 3):
    """POST with exponential backoff. Returns parsed JSON, or None on failure."""
    last_error = None
    for attempt in range(retries):
        try:
            response = requests.post(url, json=payload, params=params, timeout=timeout)
            response.raise_for_status()
            return response.json()
        except Exception as e:  # network layer: any failure is retried
            last_error = e
            time.sleep(1.5 * (attempt + 1))
    print(f"[synthetic_lethality] POST failed after {retries} attempts ({url}): {last_error}")
    return None


def _copy_number_from_local_file(gene: str, data_lake_path: str, mode: str) -> dict | None:
    """Read amplification / deep-deletion status from a local DepMap OmicsCNGene.csv if present."""
    import pandas as pd

    for filename in ("DepMap_OmicsCNGene.csv", "OmicsCNGene.csv"):
        path = os.path.join(data_lake_path, filename)
        if not os.path.exists(path):
            continue

        header = pd.read_csv(path, nrows=0)
        id_column = header.columns[0]
        column = {c.split(" (")[0].strip().upper(): c for c in header.columns}.get(gene.upper())
        if column is None:
            return None
        values = pd.read_csv(path, usecols=[id_column, column], index_col=0)[column].dropna()

        # The file is expected in log2(relative CN + 1) space (neutral = 1.0). If the column is not
        # centred there the thresholds below would be meaningless, so say so instead of guessing.
        median = float(values.median())
        scale_note = ""
        if not 0.5 <= median <= 1.5:
            scale_note = (
                f" WARNING: median {gene.upper()} value is {median:.2f}, not ~1.0 - {filename} may not be in "
                "log2(relative CN + 1) space and the amplification/deletion thresholds may not apply"
            )

        if mode == "amplification":
            altered = values[values >= DEPMAP_CN_AMPLIFIED]
            label = "AMP"
        else:
            altered = values[values <= DEPMAP_CN_DELETED]
            label = "DEL"
        return {
            "altered": {model: f"{label} (log2CN={v:.2f})" for model, v in altered.items()},
            "profiled_models": set(values.index),
            "key": "ModelID",
            "source": (
                f"local {filename} ({_file_provenance(path)}); {len(values)} models, "
                f"threshold {'>=' if mode == 'amplification' else '<='} "
                f"{DEPMAP_CN_AMPLIFIED if mode == 'amplification' else DEPMAP_CN_DELETED}{scale_note}"
            ),
        }
    return None


def _copy_number_from_cbioportal(gene: str, mode: str) -> dict | None:
    """Fetch discrete copy-number calls for ``gene`` from the cBioPortal CCLE study.

    Samples returned with alteration 0 are profiled-and-neutral, which is what makes a proper
    "neutral" control group possible rather than assuming absence of a call means no alteration.
    """
    cache_key = (gene.upper(), mode)
    if cache_key in _CN_CACHE:
        return _CN_CACHE[cache_key]

    gene_info = _http_get(f"{CBIOPORTAL_API}/genes/{urllib.parse.quote(gene.upper())}")
    if not gene_info or "entrezGeneId" not in gene_info:
        return None

    records = _http_post(
        f"{CBIOPORTAL_API}/molecular-profiles/{CCLE_STUDY_ID}_cna/discrete-copy-number/fetch",
        payload={"entrezGeneIds": [gene_info["entrezGeneId"]], "sampleListId": f"{CCLE_STUDY_ID}_all"},
        params={"discreteCopyNumberEventType": "ALL", "projection": "SUMMARY"},
    )
    if records is None:
        return None

    target = CN_AMPLIFICATION if mode == "amplification" else CN_DEEP_DELETION
    label = "AMP" if mode == "amplification" else "DEL"
    altered: dict[str, str] = {}
    profiled: set[str] = set()
    for record in records:
        key = _normalize_cell_line_name(str(record.get("sampleId", "")).split("_")[0])
        profiled.add(key)
        if record.get("alteration") == target:
            altered[key] = label

    result = {
        "altered": altered,
        "profiled_models": profiled,
        "key": "normalized_name",
        "source": (
            f"cBioPortal {CCLE_STUDY_ID} discrete CNA calls ({len(altered)} {gene.upper()} "
            f"{'amplifications' if mode == 'amplification' else 'deep deletions'} across {len(profiled)} "
            "profiled cell lines)"
        ),
    }
    _CN_CACHE[cache_key] = result
    return result


def _annotate_copy_number_status(
    models, gene: str, data_lake_path: str, mode: str, genotype_csv_path: str | None = None
):
    """Attach MUT (= altered) / WT (= neutral) / UNKNOWN copy-number status for ``gene``.

    Group labels reuse the mutation column names so that everything downstream - group selection,
    testing, reporting - is identical whatever the stratifying event is.
    Precedence: user-supplied CSV > local DepMap OmicsCNGene file > cBioPortal CCLE API.
    """
    import pandas as pd

    models = models.copy()
    models["_norm_name"] = models["StrippedCellLineName"].map(_normalize_cell_line_name)

    annotation = None
    if genotype_csv_path:
        table = pd.read_csv(genotype_csv_path)
        required = {"ModelID", "HugoSymbol"}
        if not required.issubset(table.columns):
            raise ValueError(f"{genotype_csv_path} must contain columns {sorted(required)}")
        hits = table[table["HugoSymbol"].astype(str).str.upper() == gene.upper()]
        if "Alteration" in table.columns:
            keyword = "AMP" if mode == "amplification" else "DEL"
            hits = hits[hits["Alteration"].astype(str).str.upper().str.contains(keyword, na=False)]
        annotation = {
            "altered": {r["ModelID"]: str(r.get("Alteration", mode)) for _, r in hits.iterrows()},
            "profiled_models": set(table["ModelID"].unique()),
            "key": "ModelID",
            "source": f"user-supplied genotype table {genotype_csv_path}",
        }
    if annotation is None:
        annotation = _copy_number_from_local_file(gene, data_lake_path, mode)
    if annotation is None:
        annotation = _copy_number_from_cbioportal(gene, mode)
    if annotation is None:
        raise RuntimeError(
            f"Could not determine {gene} copy-number status: no local OmicsCNGene.csv and the cBioPortal "
            "API is unreachable. Supply `mutation_csv_path` with columns ModelID,HugoSymbol[,Alteration]."
        )

    key_column = "ModelID" if annotation["key"] == "ModelID" else "_norm_name"

    def classify(key):
        if key in annotation["altered"]:
            return "MUT"
        return "WT" if key in annotation["profiled_models"] else "UNKNOWN"

    models["MutationStatus"] = models[key_column].map(classify)
    models["Variant"] = models[key_column].map(lambda k: annotation["altered"].get(k, ""))
    return models, annotation["source"]


# stratify_by value -> (altered-group label, control-group label, event name, short labels)
_STRATIFY_MODES = {
    "mutation": ("MUTANT", "WILD-TYPE", "mutation", ("MUT", "WT")),
    "amplification": ("AMPLIFIED", "NEUTRAL", "amplification", ("AMP", "NEUT")),
    "deletion": ("DEEP-DELETED", "NEUTRAL", "deep deletion", ("DEL", "NEUT")),
}
_STRATIFY_ALIASES = {
    "mut": "mutation",
    "mutation": "mutation",
    "mutated": "mutation",
    "amp": "amplification",
    "amplification": "amplification",
    "amplified": "amplification",
    "del": "deletion",
    "deletion": "deletion",
    "deleted": "deletion",
    "deep deletion": "deletion",
    "loss": "deletion",
}


def _benjamini_hochberg(pvalues):
    """Return BH-FDR adjusted q-values for a 1-D array of p-values."""
    import numpy as np

    pvalues = np.asarray(pvalues, dtype=float)
    n = len(pvalues)
    order = np.argsort(pvalues)
    ranked = pvalues[order]
    qvalues = ranked * n / (np.arange(n) + 1)
    qvalues = np.minimum.accumulate(qvalues[::-1])[::-1]
    out = np.empty(n, dtype=float)
    out[order] = np.clip(qvalues, 0, 1)
    return out


def _parse_gene_list(genes) -> list[str]:
    """Accept a list, a comma separated string, or a newline separated string of gene symbols."""
    if genes is None:
        return []
    if isinstance(genes, str):
        parts = re.split(r"[,\s]+", genes)
    else:
        parts = list(genes)
    return [str(g).strip().upper() for g in parts if str(g).strip()]


# ---------------------------------------------------------------------------
# Tool 1: data-driven synthetic lethality discovery
# ---------------------------------------------------------------------------
def discover_synthetic_lethal_candidates(
    cancer_type: str,
    target_mutation: str,
    data_lake_path: str | None = None,
    mutation_csv_path: str | None = None,
    top_n: int = 20,
    p_threshold: float = 0.05,
    fdr_threshold: float = 0.25,
    min_effect_difference: float = -0.2,
    max_mutant_mean_effect: float = -0.3,
    exclude_pan_essential: bool = True,
    output_csv_path: str | None = None,
) -> str:
    """Discover candidate synthetic-lethal partner genes of a driver mutation in a cancer type.

    Cell lines of the requested cancer type are split into mutant and wild-type groups for
    ``target_mutation``, and every gene in the DepMap CRISPR knockout screen is tested with a
    Welch t-test for stronger dependency (more negative Chronos gene effect) in the mutant group.
    Pan-essential genes, weak dependencies and low-powered comparisons are flagged or removed so
    that the surviving candidates represent genotype-selective vulnerabilities rather than
    generic essential genes.

    Parameters
    ----------
    cancer_type : str
        Cancer context to analyse, e.g. "Pancreatic Cancer", "Lung", "Colorectal Adenocarcinoma".
        Use "pan-cancer" to analyse all lineages together.
    target_mutation : str
        HUGO symbol of the mutated driver gene defining the two groups, e.g. "KRAS", "TP53", "BRCA1".
    data_lake_path : str, optional
        Directory holding DepMap_CRISPRGeneEffect.csv and DepMap_Model.csv
        (default: "./data/biomni_data/data_lake").
    mutation_csv_path : str, optional
        Path to a mutation table (columns: ModelID, HugoSymbol[, ProteinChange]) overriding the
        default mutation source. When omitted, a local DepMap somatic mutation file is used if
        present, otherwise mutation calls are fetched from the cBioPortal CCLE study.
    top_n : int, optional
        Number of top candidates to report in detail (default: 20).
    p_threshold : float, optional
        Uncorrected Welch t-test p-value cutoff (default: 0.05).
    fdr_threshold : float, optional
        Benjamini-Hochberg q-value cutoff (default: 0.25). Set to 1.0 to disable FDR filtering.
    min_effect_difference : float, optional
        Required difference (mutant mean - wild-type mean) in gene effect; must be negative
        (default: -0.2, i.e. the mutant group must be at least 0.2 more depleted).
    max_mutant_mean_effect : float, optional
        The mutant group mean gene effect must be below this value so that statistically
        significant but biologically irrelevant genes are dropped (default: -0.3).
    exclude_pan_essential : bool, optional
        Drop genes depleted in >80% of all screened cell lines (common essentials) (default: True).
    output_csv_path : str, optional
        If given, the full ranked candidate table is written to this CSV path.

    Returns
    -------
    str
        A research log containing the cohort composition, provenance, QC warnings, and a ranked
        candidate table with effect sizes, p/q-values, selectivity and pan-essentiality metrics.

    """
    import numpy as np
    import pandas as pd
    from scipy import stats

    log = [
        "=" * 78,
        f"SYNTHETIC LETHALITY DISCOVERY - {target_mutation.upper()} in {cancer_type}",
        f"Run at {datetime.now(tz=_UTC).strftime('%Y-%m-%d %H:%M UTC')}",
        "=" * 78,
        "",
    ]

    try:
        bundle = _load_depmap(data_lake_path)
    except FileNotFoundError as e:
        return f"FAILURE: {e}"

    gene_effect = bundle["gene_effect"]
    models = bundle["model"]

    # --- Step 1: cohort definition -----------------------------------------------------------
    cohort, match_note = _select_cancer_models(models, cancer_type)
    cohort = cohort[cohort["ModelID"].isin(gene_effect.index)]
    if len(cohort) == 0:
        return (
            f"FAILURE: no DepMap cell line with CRISPR data matched cancer_type='{cancer_type}'. "
            f"Available lineages: {sorted(models['OncotreeLineage'].dropna().unique().tolist())}"
        )

    log.append("STEP 1 | Cohort definition")
    log.append(f"  Cancer context: {match_note}")
    log.append(f"  Cell lines with CRISPR gene-effect data: {len(cohort)}")

    # --- Step 2: mutation stratification -----------------------------------------------------
    try:
        cohort, mutation_source = _annotate_mutation_status(
            cohort, target_mutation, bundle["data_lake_path"], mutation_csv_path
        )
    except (RuntimeError, ValueError) as e:
        return f"FAILURE: {e}"

    mutant_ids = cohort.loc[cohort["MutationStatus"] == "MUT", "ModelID"].tolist()
    wildtype_ids = cohort.loc[cohort["MutationStatus"] == "WT", "ModelID"].tolist()
    unknown_ids = cohort.loc[cohort["MutationStatus"] == "UNKNOWN", "ModelID"].tolist()

    variant_counts = cohort.loc[cohort["MutationStatus"] == "MUT", "Variant"].value_counts().head(6).to_dict()

    log.append("")
    log.append("STEP 2 | Mutation stratification")
    log.append(f"  Mutation source: {mutation_source}")
    log.append(f"  {target_mutation.upper()}-mutant lines : {len(mutant_ids)}")
    log.append(f"  {target_mutation.upper()}-wild-type lines: {len(wildtype_ids)}")
    log.append(f"  Not profiled (excluded)   : {len(unknown_ids)}")
    if variant_counts:
        log.append(f"  Variant spectrum: {variant_counts}")

    if len(mutant_ids) < 3 or len(wildtype_ids) < 3:
        log.append("")
        log.append(
            f"FAILURE: insufficient group sizes (mutant n={len(mutant_ids)}, wild-type n={len(wildtype_ids)}); "
            "at least 3 per group are required. Broaden `cancer_type` (e.g. 'pan-cancer') or supply a "
            "mutation table with wider coverage."
        )
        return "\n".join(log)

    # --- Step 3: differential dependency testing ---------------------------------------------
    mutant_matrix = gene_effect.loc[mutant_ids]
    wildtype_matrix = gene_effect.loc[wildtype_ids]
    usable_genes = mutant_matrix.columns[(mutant_matrix.notna().sum() >= 3) & (wildtype_matrix.notna().sum() >= 3)]
    mutant_matrix = mutant_matrix[usable_genes]
    wildtype_matrix = wildtype_matrix[usable_genes]

    tstat, pvalue = stats.ttest_ind(
        mutant_matrix.values, wildtype_matrix.values, axis=0, equal_var=False, nan_policy="omit"
    )
    mutant_mean = np.asarray(mutant_matrix.mean())
    wildtype_mean = np.asarray(wildtype_matrix.mean())
    pooled_sd = np.sqrt((np.asarray(mutant_matrix.std(ddof=1)) ** 2 + np.asarray(wildtype_matrix.std(ddof=1)) ** 2) / 2)
    with np.errstate(divide="ignore", invalid="ignore"):
        cohens_d = np.where(pooled_sd > 0, (mutant_mean - wildtype_mean) / pooled_sd, np.nan)

    pvalue = np.asarray(pvalue, dtype=float)
    pvalue = np.where(np.isfinite(pvalue), pvalue, 1.0)

    # Pan-essentiality across the whole screen (all lineages), used as a falsification filter.
    depleted_fraction = (gene_effect[usable_genes] < DEPLETION_THRESHOLD).sum() / gene_effect[
        usable_genes
    ].notna().sum()
    # Selectivity inside the cohort: fraction of mutant lines that are actually dependent.
    mutant_dependent_fraction = (mutant_matrix < DEPLETION_THRESHOLD).sum() / mutant_matrix.notna().sum()
    wildtype_dependent_fraction = (wildtype_matrix < DEPLETION_THRESHOLD).sum() / wildtype_matrix.notna().sum()

    results = pd.DataFrame(
        {
            "gene": [c.split(" (")[0] for c in usable_genes],
            "mutant_mean_effect": mutant_mean,
            "wildtype_mean_effect": wildtype_mean,
            "effect_difference": mutant_mean - wildtype_mean,
            "cohens_d": cohens_d,
            "t_statistic": np.asarray(tstat, dtype=float),
            "p_value": pvalue,
            "q_value": _benjamini_hochberg(pvalue),
            "pct_mutant_dependent": np.asarray(mutant_dependent_fraction) * 100,
            "pct_wildtype_dependent": np.asarray(wildtype_dependent_fraction) * 100,
            "pct_all_lines_dependent": np.asarray(depleted_fraction) * 100,
        }
    )
    results["pan_essential"] = results["pct_all_lines_dependent"] >= PAN_ESSENTIAL_FRACTION * 100

    log.append("")
    log.append("STEP 3 | Differential dependency testing (Welch t-test, mutant vs wild-type)")
    log.append(f"  Genes tested: {len(results)}")
    log.append(f"  Nominally significant (p < {p_threshold}): {(results['p_value'] < p_threshold).sum()}")
    log.append(f"  BH q < {fdr_threshold}: {(results['q_value'] < fdr_threshold).sum()}")

    # --- Step 4: filtering / QC ---------------------------------------------------------------
    selected = results[
        (results["p_value"] < p_threshold)
        & (results["q_value"] < fdr_threshold)
        & (results["effect_difference"] <= min_effect_difference)
        & (results["mutant_mean_effect"] <= max_mutant_mean_effect)
    ].copy()
    n_before_pan_essential = len(selected)
    if exclude_pan_essential:
        selected = selected[~selected["pan_essential"]]

    selected = selected.sort_values(["effect_difference", "p_value"]).reset_index(drop=True)
    selected.insert(0, "rank", np.arange(1, len(selected) + 1))

    log.append("")
    log.append("STEP 4 | Filtering and quality control")
    log.append(f"  After effect-size + significance filters: {n_before_pan_essential}")
    log.append(
        f"  Pan-essential genes removed (depleted in >{PAN_ESSENTIAL_FRACTION:.0%} of all "
        f"{gene_effect.shape[0]} lines): {n_before_pan_essential - len(selected)}"
    )
    log.append(f"  Final candidate count: {len(selected)}")

    # Positive control: the driver itself should be a selective dependency in its own mutant lines.
    control_rows = results[results["gene"].str.upper() == target_mutation.upper()]
    if len(control_rows):
        control = control_rows.iloc[0]
        passed = (control["effect_difference"] <= min_effect_difference) and (control["p_value"] < p_threshold)
        in_final = target_mutation.upper() in set(selected["gene"].str.upper())
        status = "PASS" if passed else "FAIL"
        log.append(
            f"  Positive control ({status}): {target_mutation.upper()} itself shows "
            f"delta={control['effect_difference']:.3f}, p={control['p_value']:.2e}, q={control['q_value']:.3f} "
            f"(mutant {control['mutant_mean_effect']:.3f} vs wild type {control['wildtype_mean_effect']:.3f})."
        )
        if passed and not in_final:
            log.append(
                "    Oncogene addiction is reproduced, but the control does not survive the FDR/effect filters - "
                "expected with a small wild-type group, and a reminder that the q-value cutoff is conservative here."
            )
        elif not passed:
            log.append(
                f"    WARNING: {target_mutation.upper()} is not recovered as a selective dependency of its own "
                "mutant lines. The mutation calls or the cohort may be mis-specified; treat all candidates below "
                "as unreliable until this is resolved."
            )
    else:
        log.append(f"  Positive control unavailable: {target_mutation.upper()} is not in the CRISPR library.")

    warnings = []
    if len(wildtype_ids) < 8:
        warnings.append(
            f"Low statistical power: only {len(wildtype_ids)} wild-type lines. p-values are unstable and a "
            "single outlier line can drive a candidate; replicate in a pan-cancer or external cohort."
        )
    if len(unknown_ids):
        warnings.append(
            f"{len(unknown_ids)} lines lacked mutation calls and were excluded rather than assumed wild type."
        )
    lineages = cohort["OncotreeLineage"].dropna().nunique()
    if lineages > 1:
        warnings.append(
            f"The cohort spans {lineages} lineages; lineage composition may confound the mutant/WT contrast."
        )
    warnings.append(
        "Single-perturbation correlative evidence only - no double-knockout or isogenic validation is "
        "included, so these are hypotheses, not confirmed SL interactions."
    )

    log.append("")
    log.append("QC WARNINGS (falsification checklist)")
    for warning in warnings:
        log.append(f"  - {warning}")

    # --- Step 5: report -----------------------------------------------------------------------
    log.append("")
    log.append(f"TOP {min(top_n, len(selected))} CANDIDATE SYNTHETIC LETHAL PARTNERS")
    log.append(
        f"{'rank':<5}{'gene':<12}{'MUT mean':>10}{'WT mean':>10}{'delta':>9}"
        f"{'d':>7}{'p':>11}{'q':>9}{'%MUT dep':>10}{'%all dep':>10}"
    )
    log.append("-" * 93)
    for _, row in selected.head(top_n).iterrows():
        log.append(
            f"{int(row['rank']):<5}{row['gene']:<12}{row['mutant_mean_effect']:>10.3f}"
            f"{row['wildtype_mean_effect']:>10.3f}{row['effect_difference']:>9.3f}"
            f"{row['cohens_d']:>7.2f}{row['p_value']:>11.2e}{row['q_value']:>9.3f}"
            f"{row['pct_mutant_dependent']:>9.0f}%{row['pct_all_lines_dependent']:>9.0f}%"
        )

    candidate_genes = [g for g in selected["gene"].head(top_n).tolist() if g.upper() != target_mutation.upper()]
    log.append("")
    log.append(f"CANDIDATE_GENES: {', '.join(candidate_genes)}")
    log.append(
        "  (pass this list to validate_sl_candidates_with_pubmed and analyze_ppi_network_for_sl "
        "for literature and protein-network evidence)"
    )

    if output_csv_path:
        os.makedirs(os.path.dirname(os.path.abspath(output_csv_path)), exist_ok=True)
        selected.to_csv(output_csv_path, index=False)
        log.append(f"  Full ranked table written to {output_csv_path} ({len(selected)} rows)")

    log.append("")
    log.append("PROVENANCE")
    for item in bundle["provenance"]:
        log.append(f"  - {item}")
    log.append(f"  - Mutation calls: {mutation_source}")
    log.append(
        f"  - Statistics: Welch two-sample t-test on Chronos gene effect, BH-FDR across {len(results)} genes; "
        f"depletion threshold {DEPLETION_THRESHOLD}, pan-essential cutoff {PAN_ESSENTIAL_FRACTION:.0%}"
    )
    log.append(
        f"  - Mutant models: {', '.join(cohort.loc[cohort['MutationStatus'] == 'MUT', 'StrippedCellLineName'][:40])}"
    )
    log.append(
        f"  - Wild-type models: {', '.join(cohort.loc[cohort['MutationStatus'] == 'WT', 'StrippedCellLineName'])}"
    )
    return "\n".join(log)


# ---------------------------------------------------------------------------
# Tool 2: PubMed literature validation
# ---------------------------------------------------------------------------
_SUPPORT_TERMS = [
    "synthetic lethal",
    "synthetic lethality",
    "synthetically lethal",
    "selective dependency",
    "collateral vulnerability",
    "sensitizes",
    "sensitized",
    "sensitivity to",
    "vulnerability",
    "co-dependency",
    "codependency",
    "essential in",
    "required for the survival",
    "growth inhibition",
    "impaired proliferation",
]
_REFUTATION_TERMS = [
    "not synthetic lethal",
    "no synthetic lethal",
    "failed to replicate",
    "failed to reproduce",
    "could not be reproduced",
    "not reproducible",
    "irreproducible",
    "no evidence",
    "did not confirm",
    "was not confirmed",
    "off-target",
    "independent of kras",
    "is dispensable",
    "dispensable for",
    "no significant difference",
    "contrary to",
    "challenge the",
    "call into question",
    "not required for",
    "nonessential",
    "non-essential",
    "is not essential",
    "were not essential",
    "unable to confirm",
    "could not confirm",
    "reevaluation",
    "re-evaluation",
    "reassessment",
    "controversial",
    "questioned",
    "does not depend",
    "no correlation between",
]
_DIRECT_EVIDENCE_TERMS = [
    "crispr screen",
    "shrna screen",
    "rnai screen",
    "isogenic",
    "knockout",
    "knockdown",
    "xenograft",
    "in vivo",
    "combination treatment",
    "double knockout",
]


def _entrez_esearch(query: str, max_papers: int, email: str | None, api_key: str | None) -> list[str]:
    params = {
        "db": "pubmed",
        "term": query,
        "retmax": str(max_papers),
        "retmode": "json",
        "sort": "relevance",
        "tool": "biomni_synthetic_lethality",
    }
    if email:
        params["email"] = email
    if api_key:
        params["api_key"] = api_key
    payload = _http_get(f"{ENTREZ_API}/esearch.fcgi", params=params, timeout=45)
    if not payload:
        return []
    return payload.get("esearchresult", {}).get("idlist", [])


def _entrez_efetch(pmids: list[str], email: str | None, api_key: str | None) -> list[dict]:
    """Fetch PubMed records as XML and parse title / year / journal / abstract."""
    if not pmids:
        return []
    params = {
        "db": "pubmed",
        "id": ",".join(pmids),
        "retmode": "xml",
        "tool": "biomni_synthetic_lethality",
    }
    if email:
        params["email"] = email
    if api_key:
        params["api_key"] = api_key
    xml_text = _http_get(f"{ENTREZ_API}/efetch.fcgi", params=params, timeout=60, as_json=False)
    if not xml_text:
        return []

    try:
        root = ET.fromstring(xml_text)
    except ET.ParseError as e:
        print(f"[synthetic_lethality] failed to parse PubMed XML: {e}")
        return []

    records = []
    for article in root.findall(".//PubmedArticle"):
        pmid_node = article.find(".//PMID")
        title_node = article.find(".//ArticleTitle")

        # Structured abstracts are split over several AbstractText nodes with a Label attribute.
        abstract_parts = []
        for node in article.findall(".//Abstract/AbstractText"):
            text = "".join(node.itertext()).strip()
            if not text:
                continue
            label = node.get("Label")
            abstract_parts.append(f"{label}: {text}" if label else text)

        year = None
        for path in (".//ArticleDate/Year", ".//PubDate/Year", ".//PubMedPubDate/Year"):
            node = article.find(path)
            if node is not None and node.text:
                year = node.text
                break
        if year is None:
            medline_date = article.find(".//PubDate/MedlineDate")
            if medline_date is not None and medline_date.text:
                match = re.search(r"\d{4}", medline_date.text)
                year = match.group(0) if match else None

        journal_node = article.find(".//Journal/Title")
        pub_types = [t.text for t in article.findall(".//PublicationType") if t.text]

        records.append(
            {
                "pmid": pmid_node.text if pmid_node is not None else "",
                "title": "".join(title_node.itertext()).strip() if title_node is not None else "",
                "abstract": " ".join(abstract_parts),
                "year": int(year) if year and year.isdigit() else None,
                "journal": journal_node.text if journal_node is not None else "",
                "publication_types": pub_types,
            }
        )
    return records


def _score_literature(records: list[dict], candidate_gene: str, mutated_gene: str, disease: str) -> dict:
    """Score literature support for one candidate on a 0-100 scale with an explicit breakdown."""
    support_hits, refute_hits, direct_hits = [], [], []
    disease_hits = 0
    cooccurrence = 0
    recent = 0
    current_year = datetime.now(tz=_UTC).year
    disease_tokens = [t for t in re.split(r"[^a-z]+", disease.lower()) if len(t) > 4]

    for record in records:
        text = f"{record['title']} {record['abstract']}".lower()
        if candidate_gene.lower() in text and mutated_gene.lower() in text:
            cooccurrence += 1
        if disease_tokens and any(token[:7] in text for token in disease_tokens):
            disease_hits += 1
        if record["year"] and record["year"] >= current_year - 5:
            recent += 1
        for term in _SUPPORT_TERMS:
            if term in text:
                support_hits.append((record["pmid"], term))
        for term in _REFUTATION_TERMS:
            if term in text:
                refute_hits.append((record["pmid"], term))
        for term in _DIRECT_EVIDENCE_TERMS:
            if term in text:
                direct_hits.append((record["pmid"], term))

    n = len(records)
    volume_score = min(25, n * 5)
    cooccurrence_score = min(25, cooccurrence * 8)
    support_score = min(20, len({p for p, _ in support_hits}) * 7)
    direct_score = min(15, len({p for p, _ in direct_hits}) * 5)
    context_score = min(10, disease_hits * 4)
    recency_score = min(5, recent * 2)
    refutation_penalty = min(30, len({p for p, _ in refute_hits}) * 10)

    total = max(
        0,
        volume_score
        + cooccurrence_score
        + support_score
        + direct_score
        + context_score
        + recency_score
        - refutation_penalty,
    )

    if n == 0:
        interpretation = "UNEXPLORED - no PubMed record links this pair; novel but entirely unvalidated."
    elif refutation_penalty >= 20:
        interpretation = "CONTESTED - the literature contains explicit negative or non-replication statements."
    elif refutation_penalty > 0 and total >= 45:
        interpretation = (
            "SUPPORTED BUT CONTESTED - widely reported, yet at least one paper reports a negative or "
            "non-replication result; the interaction must be treated as unresolved, not established."
        )
    elif refutation_penalty > 0:
        interpretation = (
            "CONTESTED-SPARSE - little supporting literature and at least one negative/non-replication "
            "statement; verify the flagged PMIDs before spending experimental effort."
        )
    elif total >= 65:
        interpretation = "WELL SUPPORTED - repeatedly reported, including direct experimental evidence."
    elif total >= 35:
        interpretation = "EMERGING - some supporting reports, but sparse or indirect."
    else:
        interpretation = "WEAK - mentioned in the literature without specific support for this interaction."

    return {
        "n_papers": n,
        "score": total,
        "interpretation": interpretation,
        "breakdown": {
            "volume": volume_score,
            "gene_pair_cooccurrence": cooccurrence_score,
            "sl_language": support_score,
            "direct_experimental": direct_score,
            "disease_context": context_score,
            "recency": recency_score,
            "refutation_penalty": -refutation_penalty,
        },
        "support_terms": sorted({t for _, t in support_hits}),
        "refutation_terms": sorted({t for _, t in refute_hits}),
        "direct_evidence_terms": sorted({t for _, t in direct_hits}),
        "refuting_pmids": sorted({p for p, _ in refute_hits}),
    }


def validate_sl_candidates_with_pubmed(
    disease: str,
    mutated_gene: str,
    candidate_genes: list[str] | str,
    max_papers_per_gene: int = 8,
    min_year: int | None = None,
    include_abstracts: bool = True,
    abstract_chars: int = 700,
    email: str | None = None,
    api_key: str | None = None,
) -> str:
    """Validate synthetic-lethal candidate genes against the PubMed literature via the NCBI Entrez API.

    For every candidate gene a PubMed query combining the disease, the mutated driver gene and the
    candidate is executed with esearch; the matching records are retrieved with efetch and their
    abstracts are parsed from XML. Each candidate receives a 0-100 literature support score whose
    components are reported explicitly, including a penalty for refuting/non-replication language,
    so that plausible-but-fragile candidates are actively down-ranked.

    Parameters
    ----------
    disease : str
        Disease context used in the query, e.g. "Pancreatic cancer".
    mutated_gene : str
        The mutated driver gene, e.g. "KRAS".
    candidate_genes : list[str] | str
        Candidate SL partner genes; a list, or a comma/whitespace separated string.
    max_papers_per_gene : int, optional
        Maximum number of PubMed records to retrieve per candidate (default: 8).
    min_year : int, optional
        Only consider papers published in or after this year.
    include_abstracts : bool, optional
        Include the (truncated) abstract text of the top paper per candidate (default: True).
    abstract_chars : int, optional
        Number of abstract characters to show when ``include_abstracts`` is True (default: 700).
    email : str, optional
        Contact e-mail passed to Entrez, as recommended by NCBI usage policy.
    api_key : str, optional
        NCBI API key; raises the rate limit from 3 to 10 requests/second.

    Returns
    -------
    str
        A research log with, per candidate, the PubMed query used, the number of hits, the support
        score with its breakdown, refuting evidence, and the key papers (PMID, year, journal, title,
        optional abstract), followed by a ranked summary table.

    """
    genes = _parse_gene_list(candidate_genes)
    if not genes:
        return "FAILURE: no candidate genes supplied."

    api_key = api_key or os.getenv("NCBI_API_KEY")
    email = email or os.getenv("NCBI_EMAIL")
    delay = 0.12 if api_key else 0.4  # NCBI: 10 req/s with a key, 3 req/s without

    log = [
        "=" * 78,
        f"PUBMED LITERATURE VALIDATION - {mutated_gene.upper()} SL partners in {disease}",
        f"Candidates: {', '.join(genes)}",
        f"Run at {datetime.now(tz=_UTC).strftime('%Y-%m-%d %H:%M UTC')}",
        "=" * 78,
    ]

    date_filter = f' AND ("{min_year}"[PDAT] : "3000"[PDAT])' if min_year else ""
    summaries = []

    for gene in genes:
        pair = f'"{mutated_gene}"[Title/Abstract] AND "{gene}"[Title/Abstract]{date_filter}'
        # Tier 1 anchors the disease context; tier 2 deliberately hunts for synthetic-lethality
        # claims AND their refutations outside that context; tier 3 is the unrestricted pair.
        tiers = [
            (
                "disease-specific",
                f'("{disease}"[Title/Abstract] OR "{disease}"[MeSH Terms]) AND {pair}',
            ),
            (
                "SL-focused",
                f"{pair} AND (synthetic lethal*[Title/Abstract] OR dependency[Title/Abstract] OR "
                "vulnerability[Title/Abstract] OR sensitiz*[Title/Abstract] OR screen[Title/Abstract])",
            ),
            (
                "refutation-focused",
                f"{pair} AND (nonessential[Title/Abstract] OR non-essential[Title/Abstract] OR "
                "dispensable[Title/Abstract] OR reproducib*[Title/Abstract] OR off-target[Title/Abstract] OR "
                "controvers*[Title/Abstract] OR reevaluat*[Title/Abstract] OR re-evaluat*[Title/Abstract] OR "
                "reassess*[Title/Abstract] OR fail*[Title/Abstract])",
            ),
        ]

        pmids: list[str] = []
        used_queries = []
        for label, query in tiers:
            hits = _entrez_esearch(query, max_papers_per_gene, email, api_key)
            time.sleep(delay)
            new = [p for p in hits if p not in pmids]
            used_queries.append(f"[{label}] {len(hits)} hit(s), {len(new)} new")
            pmids.extend(new)

        broadened = False
        if not pmids:
            # SL interactions are often first reported outside the disease context of interest.
            pmids = _entrez_esearch(pair, max_papers_per_gene, email, api_key)
            used_queries.append(f"[gene-pair fallback] {len(pmids)} hit(s)")
            broadened = True
            time.sleep(delay)

        pmids = pmids[: max_papers_per_gene * 2]
        records = _entrez_efetch(pmids, email, api_key)
        time.sleep(delay)
        scored = _score_literature(records, gene, mutated_gene, disease)
        scored["gene"] = gene
        scored["query"] = pair
        scored["broadened"] = broadened
        summaries.append(scored)

        log.append("")
        log.append(f"### {gene}")
        log.append(f"  Queries: {'; '.join(used_queries)}")
        log.append(f"  Base pair term: {pair}")
        if broadened:
            log.append("  NOTE: no disease-specific or SL-focused paper found; the unrestricted pair query was used.")
        log.append(f"  Records retrieved: {scored['n_papers']}")
        log.append(f"  Literature support score: {scored['score']}/100 -> {scored['interpretation']}")
        log.append(f"  Score breakdown: {json.dumps(scored['breakdown'])}")
        if scored["support_terms"]:
            log.append(f"  Supporting language: {', '.join(scored['support_terms'][:8])}")
        if scored["direct_evidence_terms"]:
            log.append(f"  Direct experimental evidence terms: {', '.join(scored['direct_evidence_terms'][:8])}")
        if scored["refutation_terms"]:
            log.append(
                f"  REFUTING language: {', '.join(scored['refutation_terms'][:8])} "
                f"(PMIDs {', '.join(scored['refuting_pmids'][:5])})"
            )
        else:
            log.append("  REFUTING language: none detected in the retrieved abstracts.")

        for record in records[:3]:
            log.append(
                f"  - PMID {record['pmid']} ({record['year'] or 'n.d.'}, {record['journal']}): {record['title']}"
            )
            if include_abstracts and record["abstract"]:
                abstract = record["abstract"][:abstract_chars]
                suffix = "..." if len(record["abstract"]) > abstract_chars else ""
                log.append(f"      Abstract: {abstract}{suffix}")

    summaries.sort(key=lambda s: s["score"], reverse=True)
    log.append("")
    log.append("LITERATURE SUPPORT RANKING")
    log.append(f"{'gene':<12}{'papers':>8}{'score':>8}  interpretation")
    log.append("-" * 78)
    for summary in summaries:
        log.append(f"{summary['gene']:<12}{summary['n_papers']:>8}{summary['score']:>8}  {summary['interpretation']}")

    log.append("")
    log.append("PROVENANCE")
    log.append(f"  - NCBI Entrez E-utilities ({ENTREZ_API}), db=pubmed, esearch + efetch (XML abstracts)")
    log.append(f"  - API key in use: {'yes' if api_key else 'no (3 requests/second limit)'}")
    log.append(
        "  - Scores are computed by deterministic keyword rules over retrieved abstracts, not by an LLM; "
        "a low score means 'not documented', not 'disproven'."
    )
    log.append(
        "  - Keyword-based refutation detection can fire on a phrase used in an unrelated context; always open "
        "the flagged PMIDs before acting on a CONTESTED verdict."
    )
    return "\n".join(log)


# ---------------------------------------------------------------------------
# Tool 3: STRING protein-protein interaction analysis
# ---------------------------------------------------------------------------
def analyze_ppi_network_for_sl(
    target_gene: str,
    candidate_genes: list[str] | str,
    species: int = 9606,
    required_score: int = 400,
    analyze_shared_partners: bool = True,
    partner_limit: int = 250,
    run_enrichment: bool = True,
) -> str:
    """Analyse STRING-DB protein interactions between a mutated driver gene and SL candidates.

    Physical/functional proximity is evaluated in two ways: (i) the direct STRING edge between the
    driver and each candidate, broken down by evidence channel (experimental, database, co-expression,
    text-mining), and (ii) the number of shared first-degree interaction partners, which captures
    pathway-level ("hop-2") proximity typical of genuine SL pairs that never touch physically.

    Parameters
    ----------
    target_gene : str
        The mutated driver gene, e.g. "KRAS".
    candidate_genes : list[str] | str
        Candidate SL partner genes; a list, or a comma/whitespace separated string.
    species : int, optional
        NCBI taxonomy identifier (default: 9606, Homo sapiens).
    required_score : int, optional
        Minimum STRING combined score, 0-1000 (default: 400 = medium confidence).
    analyze_shared_partners : bool, optional
        Compute shared first-degree neighbours between the driver and each candidate (default: True).
    partner_limit : int, optional
        Maximum number of interaction partners fetched per gene for the shared-neighbour analysis
        (default: 250).
    run_enrichment : bool, optional
        Run STRING functional enrichment over the driver + candidate set (default: True).

    Returns
    -------
    str
        A research log with an identifier mapping report, a direct-interaction table with per-channel
        evidence scores, shared-partner counts, a functional-proximity classification per candidate,
        and the enriched pathways of the joint network.

    """
    genes = _parse_gene_list(candidate_genes)
    if not genes:
        return "FAILURE: no candidate genes supplied."
    target_gene = target_gene.upper()
    all_genes = [target_gene] + [g for g in genes if g != target_gene]

    log = [
        "=" * 78,
        f"STRING PPI NETWORK ANALYSIS - {target_gene} vs {len(genes)} candidate(s)",
        f"Species {species}, minimum combined score {required_score}",
        f"Run at {datetime.now(tz=_UTC).strftime('%Y-%m-%d %H:%M UTC')}",
        "=" * 78,
        "",
    ]

    # --- Step 1: identifier mapping -----------------------------------------------------------
    mapping = _http_get(
        f"{STRING_API}/json/get_string_ids",
        params={"identifiers": "\r".join(all_genes), "species": species, "limit": 1, "caller_identity": "biomni"},
        timeout=60,
    )
    resolved = {}
    if mapping:
        for item in mapping:
            resolved[str(item.get("queryItem", "")).upper()] = item.get("preferredName", "")
    unresolved = [g for g in all_genes if g not in resolved]

    log.append("STEP 1 | Identifier mapping")
    log.append(f"  Resolved {len(resolved)}/{len(all_genes)} symbols in STRING v12.")
    if unresolved:
        log.append(f"  WARNING: unresolved symbols (excluded): {', '.join(unresolved)}")
    if target_gene not in resolved:
        return "\n".join(log + ["", f"FAILURE: the driver gene {target_gene} could not be mapped in STRING."])

    # --- Step 2: direct interactions ----------------------------------------------------------
    network = _http_get(
        f"{STRING_API}/json/network",
        params={
            "identifiers": "\r".join(all_genes),
            "species": species,
            "required_score": required_score,
            "caller_identity": "biomni",
        },
        timeout=90,
    )
    if network is None:
        return "\n".join(log + ["", "FAILURE: the STRING network endpoint is unreachable."])

    edges = {}
    for edge in network:
        a, b = str(edge.get("preferredName_A", "")).upper(), str(edge.get("preferredName_B", "")).upper()
        edges[frozenset((a, b))] = edge

    log.append("")
    log.append("STEP 2 | Direct interaction with the driver gene")
    log.append(f"{'candidate':<12}{'combined':>10}{'experim.':>10}{'database':>10}{'coexpr.':>10}{'textmine':>10}")
    log.append("-" * 62)

    direct = {}
    for gene in genes:
        edge = edges.get(frozenset((target_gene, resolved.get(gene, gene))))
        if edge:
            direct[gene] = {
                "combined": float(edge.get("score", 0)),
                "experimental": float(edge.get("escore", 0)),
                "database": float(edge.get("dscore", 0)),
                "coexpression": float(edge.get("ascore", 0)),
                "textmining": float(edge.get("tscore", 0)),
            }
            d = direct[gene]
            log.append(
                f"{gene:<12}{d['combined']:>10.3f}{d['experimental']:>10.3f}"
                f"{d['database']:>10.3f}{d['coexpression']:>10.3f}{d['textmining']:>10.3f}"
            )
        else:
            direct[gene] = None
            log.append(f"{gene:<12}  (no direct edge at combined score >= {required_score})")

    # --- Step 3: shared interaction partners ---------------------------------------------------
    shared_counts = {}
    if analyze_shared_partners:

        def partners_of(gene: str) -> set:
            payload = _http_get(
                f"{STRING_API}/json/interaction_partners",
                params={
                    "identifiers": gene,
                    "species": species,
                    "required_score": required_score,
                    "limit": partner_limit,
                    "caller_identity": "biomni",
                },
                timeout=90,
            )
            if not payload:
                return set()
            return {str(p.get("preferredName_B", "")).upper() for p in payload}

        target_partners = partners_of(target_gene)
        log.append("")
        log.append("STEP 3 | Shared first-degree interaction partners (pathway-level proximity)")
        log.append(f"  {target_gene} has {len(target_partners)} partners at score >= {required_score}.")
        for gene in genes:
            if gene not in resolved:
                continue
            candidate_partners = partners_of(gene)
            shared = sorted(target_partners & candidate_partners)
            shared_counts[gene] = shared
            preview = ", ".join(shared[:10]) + ("..." if len(shared) > 10 else "")
            log.append(f"  {gene:<10} partners={len(candidate_partners):<4} shared={len(shared):<4} {preview}")

    # --- Step 4: proximity classification ------------------------------------------------------
    log.append("")
    log.append("STEP 4 | Functional proximity classification")
    classifications = {}
    for gene in genes:
        edge = direct.get(gene)
        shared = shared_counts.get(gene, [])
        if edge and edge["experimental"] >= 0.4:
            verdict = "DIRECT-EXPERIMENTAL: physical/experimental interaction reported with the driver."
        elif edge and edge["combined"] >= 0.7:
            verdict = "DIRECT-HIGH-CONFIDENCE: strong STRING edge, largely from databases/text mining."
        elif edge:
            verdict = "DIRECT-WEAK: an edge exists but only at medium confidence."
        elif len(shared) >= 10:
            verdict = f"PATHWAY-PROXIMAL: no direct edge, but {len(shared)} shared partners (same functional module)."
        elif len(shared) >= 3:
            verdict = f"WEAKLY-PROXIMAL: {len(shared)} shared partners."
        else:
            verdict = (
                "DISTANT: no direct edge and few shared partners - either a parallel-pathway SL "
                "(mechanistically plausible) or a statistical artefact; prioritise mechanistic follow-up."
            )
        classifications[gene] = verdict
        log.append(f"  {gene:<12}{verdict}")

    log.append("")
    log.append(
        "  Interpretation note: absence of a PPI edge does NOT refute synthetic lethality. Classic SL pairs "
        "(e.g. BRCA1-PARP1) act through parallel pathways rather than physical binding. A strong text-mining-only "
        "edge, on the other hand, may merely reflect co-citation bias."
    )

    # --- Step 5: functional enrichment ---------------------------------------------------------
    if run_enrichment:
        enrichment = _http_get(
            f"{STRING_API}/json/enrichment",
            params={"identifiers": "\r".join(all_genes), "species": species, "caller_identity": "biomni"},
            timeout=90,
        )
        log.append("")
        log.append("STEP 5 | Functional enrichment of the driver + candidate set")
        if enrichment is None:
            log.append("  Enrichment endpoint unreachable; skipped.")
        elif len(enrichment) == 0:
            log.append(
                "  No enriched term - the candidates do not share an annotated pathway with the driver. "
                "Expected for parallel-pathway SL, but it also means no pathway-level corroboration is available."
            )
        else:
            filtered = [e for e in enrichment if e.get("category") in {"Process", "KEGG", "RCTM", "WikiPathways"}]
            filtered.sort(key=lambda e: float(e.get("fdr", 1)))
            for item in filtered[:10]:
                log.append(
                    f"  [{item.get('category')}] {item.get('description')} "
                    f"(genes {item.get('number_of_genes')}, FDR {float(item.get('fdr', 1)):.2e})"
                )
            if not filtered:
                log.append("  Enriched terms exist but none in the Process/KEGG/Reactome/WikiPathways categories.")

    log.append("")
    log.append("PROVENANCE")
    log.append(f"  - STRING-DB REST API ({STRING_API}), species={species}, required_score={required_score}")
    endpoints = ["get_string_ids", "network"]
    if analyze_shared_partners:
        endpoints.append("interaction_partners")
    if run_enrichment:
        endpoints.append("enrichment")
    log.append(f"  - Endpoints used: {', '.join(endpoints)}")
    log.append(f"  - Query identifiers: {', '.join(all_genes)}")
    return "\n".join(log)


# ---------------------------------------------------------------------------
# Tool 4: confounder / falsification analysis
# ---------------------------------------------------------------------------
def check_dependency_confounders(
    target_mutation: str,
    candidate_genes: list[str] | str,
    cancer_type: str = "pan-cancer",
    data_lake_path: str | None = None,
    mutation_csv_path: str | None = None,
    top_biomarkers: int = 5,
) -> str:
    """Test whether a candidate dependency is really driven by the mutation or by a confounder.

    A gene can look synthetic-lethal with a driver mutation simply because the two groups differ in
    something else. This tool runs three falsification checks per candidate across the whole DepMap
    panel and reports whether the driver hypothesis survives:

    1. Co-dependency: correlation between the candidate's and the driver's gene-effect profiles.
       Genuine members of the driver's pathway co-vary with it (e.g. RAF1 with KRAS, r about 0.39).
    2. Expression biomarker: the genes whose expression best predicts the candidate's dependency.
       If the top biomarker is the candidate's own paralog rather than anything driver-related, the
       dependency is a paralog-loss effect, not synthetic lethality with the driver (the classic
       VPS4A/VPS4B case).
    3. Lineage effect: how concentrated the dependency is in single lineages, which is the usual
       source of a spurious mutant-versus-wild-type difference.

    Parameters
    ----------
    target_mutation : str
        The driver gene whose mutation defined the candidates, e.g. "KRAS".
    candidate_genes : list[str] | str
        Candidate genes to scrutinise; a list, or a comma/whitespace separated string.
    cancer_type : str, optional
        Restrict the lineage effect summary to this context; correlations always use all lineages
        so that the sample size is adequate (default: "pan-cancer").
    data_lake_path : str, optional
        Directory holding the DepMap files (default: "./data/biomni_data/data_lake").
    mutation_csv_path : str, optional
        Optional mutation table overriding the default mutation source.
    top_biomarkers : int, optional
        Number of top expression biomarkers to report per candidate (default: 5).

    Returns
    -------
    str
        A research log with, per candidate, the driver co-dependency, the strongest expression
        biomarkers of the dependency, the paralog check, the lineage concentration, and a verdict of
        DRIVER-CONSISTENT / CONFOUNDED / UNEXPLAINED.

    """
    import numpy as np
    import pandas as pd
    from scipy import stats

    genes = _parse_gene_list(candidate_genes)
    if not genes:
        return "FAILURE: no candidate genes supplied."

    try:
        bundle = _load_depmap(data_lake_path)
        expression_bundle = _load_expression(data_lake_path)
    except FileNotFoundError as e:
        return f"FAILURE: {e}"

    gene_effect = bundle["gene_effect"]
    effect_columns = bundle["gene_columns"]
    expression = expression_bundle["expression"]
    driver = target_mutation.upper()

    log = [
        "=" * 78,
        f"CONFOUNDER ANALYSIS - are these dependencies really about {driver}?",
        f"Candidates: {', '.join(genes)}",
        f"Run at {datetime.now(tz=_UTC).strftime('%Y-%m-%d %H:%M UTC')}",
        "=" * 78,
        "",
    ]

    if driver not in effect_columns:
        return "\n".join(log + [f"FAILURE: {driver} is not in the CRISPR library, so co-dependency cannot be tested."])

    driver_effect = gene_effect[effect_columns[driver]]

    # Reference scale: how strongly do known members of the driver's own pathway co-vary with it?
    reference = {}
    for pathway_gene in ("RAF1", "SHOC2", "MAPK1", "BRAF", "PTPN11", "SOS1"):
        if pathway_gene in effect_columns and pathway_gene != driver:
            y = gene_effect[effect_columns[pathway_gene]]
            mask = driver_effect.notna() & y.notna()
            if mask.sum() > 30:
                reference[pathway_gene] = stats.pearsonr(driver_effect[mask], y[mask])[0]
    if reference:
        best_reference = max(reference.values())
        log.append("REFERENCE SCALE | co-dependency of known MAPK-pathway genes with " + driver)
        log.append("  " + ", ".join(f"{g} r={r:+.3f}" for g, r in sorted(reference.items(), key=lambda x: -x[1])))
        log.append(f"  Strongest reference r = {best_reference:+.3f}; candidates are judged against this scale.")
    else:
        best_reference = 0.3
        log.append(f"REFERENCE SCALE | unavailable; falling back to r = {best_reference:.2f} as the bar.")

    # Lineage context for the concentration check.
    cohort, _ = _select_cancer_models(bundle["model"], cancer_type)
    cohort_ids = set(cohort["ModelID"]) & set(gene_effect.index)
    lineage_by_model = bundle["model"].set_index("ModelID")["OncotreeLineage"]

    shared_models = gene_effect.index.intersection(expression.index)
    expression_matrix = expression.loc[shared_models]
    expression_values = expression_matrix.to_numpy(dtype=float)
    expression_valid = ~np.isnan(expression_values)

    verdicts = {}
    for gene in genes:
        log.append("")
        log.append(f"### {gene}")
        if gene not in effect_columns:
            log.append("  Not present in the CRISPR library; skipped.")
            continue

        candidate_effect = gene_effect[effect_columns[gene]]

        # --- check 1: co-dependency with the driver ---
        mask = driver_effect.notna() & candidate_effect.notna()
        codependency, codependency_p = stats.pearsonr(driver_effect[mask], candidate_effect[mask])
        log.append(
            f"  1) Co-dependency with {driver}: r={codependency:+.3f} (p={codependency_p:.1e}, n={int(mask.sum())}) "
            f"vs strongest pathway reference r={best_reference:+.3f}"
        )

        # --- check 2: expression biomarkers of the dependency ---
        y = candidate_effect.loc[shared_models].to_numpy(dtype=float)
        y_valid = ~np.isnan(y)
        # Pearson correlation of the dependency against every expressed gene, computed from raw
        # sums so that each column uses only the cell lines where both values are present.
        usable = expression_valid & y_valid[:, None]
        counts = usable.sum(axis=0).astype(float)
        x = np.where(usable, expression_values, 0.0)
        y_masked = np.where(usable, y[:, None], 0.0)
        sum_x = x.sum(axis=0)
        sum_y = y_masked.sum(axis=0)
        sum_xy = (x * y_masked).sum(axis=0)
        sum_xx = (x * x).sum(axis=0)
        sum_yy = (y_masked * y_masked).sum(axis=0)
        numerator = counts * sum_xy - sum_x * sum_y
        denominator = np.sqrt(
            np.clip(counts * sum_xx - sum_x**2, 0, None) * np.clip(counts * sum_yy - sum_y**2, 0, None)
        )
        with np.errstate(divide="ignore", invalid="ignore"):
            correlations = np.where(denominator > 0, numerator / denominator, np.nan)
        correlations = np.where(counts >= 100, correlations, np.nan)

        biomarker_series = pd.Series(correlations, index=[c.split(" (")[0] for c in expression_matrix.columns])
        strongest = biomarker_series.reindex(biomarker_series.abs().sort_values(ascending=False).index)
        strongest = strongest[strongest.index != gene].head(top_biomarkers)
        log.append(
            f"  2) Strongest expression biomarkers of {gene} dependency (positive r = high expression -> less depleted):"
        )
        for biomarker, value in strongest.items():
            log.append(f"       {biomarker:<12} r={value:+.3f}")

        driver_expression_r = biomarker_series.get(driver, float("nan"))
        driver_rank = (
            int((biomarker_series.abs() > abs(driver_expression_r)).sum()) + 1
            if pd.notna(driver_expression_r)
            else None
        )
        if driver_rank:
            log.append(
                f"       {driver} expression itself: r={driver_expression_r:+.3f} (rank {driver_rank} of "
                f"{int(biomarker_series.notna().sum())} genes)"
            )

        # Paralog check: scan the whole gene family explicitly, since a paralog that ranks just
        # outside the top biomarkers (VPS4B for VPS4A) still explains the dependency away.
        family_prefixes = {re.sub(r"(?<=\d)[A-Z]$", "", gene), re.sub(r"\d+$", "", gene)}
        family_prefixes = {p for p in family_prefixes if len(p) >= 3 and p != gene}
        family = [
            symbol
            for symbol in biomarker_series.index
            if symbol != gene and any(symbol.startswith(prefix) for prefix in family_prefixes)
        ]
        paralog_hit = None
        if family:
            family_series = biomarker_series.reindex(sorted(set(family))).dropna()
            family_series = family_series.reindex(family_series.abs().sort_values(ascending=False).index)
            shown = ", ".join(f"{symbol} r={value:+.3f}" for symbol, value in family_series.head(4).items())
            log.append(f"       Paralog family ({'/'.join(sorted(family_prefixes))}*): {shown}")
            if len(family_series) and abs(family_series.iloc[0]) >= 0.25:
                paralog_hit = (family_series.index[0], family_series.iloc[0])
        if paralog_hit:
            paralog_rank = int((biomarker_series.abs() > abs(paralog_hit[1])).sum()) + 1
            log.append(
                f"       PARALOG ALERT: {paralog_hit[0]} expression (r={paralog_hit[1]:+.3f}, rank {paralog_rank} of "
                f"{int(biomarker_series.notna().sum())}) predicts this dependency far better than {driver} "
                f"(rank {driver_rank}) - this looks like a paralog-loss dependency, independent of {driver} status."
            )

        # --- check 3: lineage concentration ---
        depleted = candidate_effect[candidate_effect < DEPLETION_THRESHOLD]
        lineage_counts = lineage_by_model.reindex(depleted.index).value_counts()
        total_depleted = int(lineage_counts.sum())
        if total_depleted:
            top_lineage = lineage_counts.index[0]
            top_share = lineage_counts.iloc[0] / total_depleted * 100
            log.append(
                f"  3) Lineage concentration: {total_depleted} dependent lines overall, "
                f"top lineage {top_lineage} holds {top_share:.0f}%"
            )
        else:
            top_share = 0.0
            log.append("  3) Lineage concentration: no line passes the depletion threshold.")

        in_cohort = [m for m in depleted.index if m in cohort_ids]
        log.append(f"     Within '{cancer_type}': {len(in_cohort)} dependent lines")

        # --- verdict ---
        if paralog_hit:
            verdict = (
                f"CONFOUNDED - the dependency tracks {paralog_hit[0]} expression, not {driver} status. "
                "Stratify by that biomarker before any experiment; as it stands this is not a "
                f"{driver} synthetic lethality."
            )
        elif abs(codependency) >= best_reference * 0.5 and codependency > 0:
            verdict = f"DRIVER-CONSISTENT - co-dependency with {driver} is on the scale of its known pathway members."
        elif pd.notna(driver_expression_r) and driver_rank and driver_rank <= 100:
            verdict = f"DRIVER-PLAUSIBLE - {driver} expression is among the strongest predictors of this dependency."
        elif top_share >= 50:
            verdict = (
                f"CONFOUNDED - {top_share:.0f}% of dependent lines come from one lineage; the mutant/wild-type "
                "contrast may just be lineage composition."
            )
        else:
            verdict = (
                f"UNEXPLAINED - no co-dependency with {driver}, no {driver}-related biomarker and no single "
                "lineage explains it. Could still be a genuine parallel-pathway SL, but the mechanism is open."
            )
        verdicts[gene] = verdict
        log.append(f"  VERDICT: {verdict}")

    log.append("")
    log.append("SUMMARY")
    for gene, verdict in verdicts.items():
        log.append(f"  {gene:<12}{verdict.split(' - ')[0]}")

    log.append("")
    log.append("PROVENANCE")
    for item in bundle["provenance"]:
        log.append(f"  - {item}")
    log.append(f"  - {expression_bundle['provenance']}")
    log.append(
        "  - Correlations are Pearson across all DepMap lines with both measurements; a confounded verdict is a "
        "reason to stratify the analysis, not proof that the interaction is absent."
    )
    return "\n".join(log)


# ---------------------------------------------------------------------------
# Tool 5: integrated evidence dossier
# ---------------------------------------------------------------------------
def _parse_candidate_line(discovery_report: str) -> list[str]:
    for line in discovery_report.splitlines():
        if line.startswith("CANDIDATE_GENES:"):
            return _parse_gene_list(line.split(":", 1)[1])
    return []


def _parse_discovery_table(discovery_report: str) -> dict:
    """Recover the numeric candidate table from a discovery report for integrated scoring."""
    stats = {}
    for line in discovery_report.splitlines():
        match = re.match(
            r"^\s*(\d+)\s+([A-Z0-9\-\.]+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+"
            r"(-?\d+\.\d+)\s+([\d.eE+-]+)\s+(\d+\.\d+)\s+(\d+)%\s+(\d+)%\s*$",
            line,
        )
        if match:
            stats[match.group(2)] = {
                "rank": int(match.group(1)),
                "mutant_mean_effect": float(match.group(3)),
                "wildtype_mean_effect": float(match.group(4)),
                "effect_difference": float(match.group(5)),
                "cohens_d": float(match.group(6)),
                "p_value": float(match.group(7)),
                "q_value": float(match.group(8)),
                "pct_mutant_dependent": float(match.group(9)),
                "pct_all_lines_dependent": float(match.group(10)),
            }
    return stats


def _parse_literature_scores(literature_report: str) -> dict:
    scores = {}
    current = None
    for line in literature_report.splitlines():
        header = re.match(r"^### ([A-Z0-9\-\.]+)\s*$", line)
        if header:
            current = header.group(1)
            continue
        if current:
            score_match = re.search(r"Literature support score: (\d+)/100 -> (.+)$", line)
            if score_match:
                scores[current] = {
                    "score": int(score_match.group(1)),
                    "interpretation": score_match.group(2).strip(),
                    "refuted": False,
                }
            if line.strip().startswith("REFUTING language:") and "none detected" not in line:
                scores.setdefault(current, {"score": 0, "interpretation": "", "refuted": False})
                scores[current]["refuted"] = True
    return scores


def _parse_ppi_classification(ppi_report: str) -> dict:
    verdicts = {}
    in_section = False
    for line in ppi_report.splitlines():
        if line.startswith("STEP 4 |"):
            in_section = True
            continue
        if in_section:
            if line.startswith("  Interpretation note") or line.startswith("STEP 5"):
                in_section = False
                continue
            match = re.match(r"^\s{2}([A-Z0-9\-\.]+)\s{2,}([A-Z\-]+):", line)
            if match:
                verdicts[match.group(1)] = match.group(2)
    return verdicts


def generate_sl_evidence_dossier(
    cancer_type: str,
    target_mutation: str,
    top_n: int = 5,
    data_lake_path: str | None = None,
    mutation_csv_path: str | None = None,
    max_papers_per_gene: int = 8,
    required_score: int = 400,
    email: str | None = None,
    api_key: str | None = None,
    check_confounders: bool = True,
    include_stage_logs: bool = False,
) -> str:
    """Run the full synthetic lethality pipeline and produce an integrated evidence dossier.

    Chains ``discover_synthetic_lethal_candidates`` (statistical evidence),
    ``validate_sl_candidates_with_pubmed`` (literature evidence), ``analyze_ppi_network_for_sl``
    (protein-network evidence) and ``check_dependency_confounders`` (falsification), then integrates
    the streams into a per-candidate dossier with a confidence grade, supporting and contradicting
    evidence, remaining uncertainties, a minimal validation experiment and a Go / Hold / No-go
    recommendation.

    The integrated dossier is returned FIRST and the verbose per-stage logs last, because agent
    frameworks (Biomni included) crop tool output to the first ~10,000 characters; putting the
    recommendation at the end would hide exactly the part that must not be missed.

    Because the three streams are not independent (all three ultimately derive from partially
    overlapping cell-line and literature resources), evidence is weighted by directness rather than
    simply summed: genotype-selective CRISPR effect size dominates, literature and PPI evidence
    modulate the confidence grade, and explicit refuting evidence can veto a Go recommendation.

    Parameters
    ----------
    cancer_type : str
        Cancer context, e.g. "Pancreatic Cancer".
    target_mutation : str
        Mutated driver gene defining the genotype contrast, e.g. "KRAS".
    top_n : int, optional
        Number of top statistical candidates carried into literature and PPI validation (default: 5).
    data_lake_path : str, optional
        Directory holding the DepMap files (default: "./data/biomni_data/data_lake").
    mutation_csv_path : str, optional
        Optional mutation table overriding the default mutation source.
    max_papers_per_gene : int, optional
        Maximum PubMed records per candidate (default: 8).
    required_score : int, optional
        Minimum STRING combined score (default: 400).
    email : str, optional
        Contact e-mail for the NCBI Entrez API.
    api_key : str, optional
        NCBI API key for a higher Entrez rate limit.
    check_confounders : bool, optional
        Run the confounder/falsification analysis and let a CONFOUNDED verdict force a No-go
        (default: True). Requires the DepMap expression matrix.
    include_stage_logs : bool, optional
        Append the full per-stage research logs after the dossier (default: False). Leave this off
        when calling from an agent, or the output will be cropped.

    Returns
    -------
    str
        The integrated evidence dossier first (confidence grade, supporting/contradicting evidence,
        proposed validation experiment and Go/Hold/No-go recommendation per candidate), optionally
        followed by the per-stage research logs.

    """
    sections = []

    discovery = discover_synthetic_lethal_candidates(
        cancer_type=cancer_type,
        target_mutation=target_mutation,
        data_lake_path=data_lake_path,
        mutation_csv_path=mutation_csv_path,
        top_n=max(top_n, 10),
    )
    sections.append(discovery)
    if discovery.startswith("FAILURE") or "CANDIDATE_GENES:" not in discovery:
        return "\n\n".join(["Pipeline stopped: candidate discovery produced no usable candidates.", discovery])

    candidates = _parse_candidate_line(discovery)[:top_n]
    discovery_stats = _parse_discovery_table(discovery)

    literature = validate_sl_candidates_with_pubmed(
        disease=cancer_type,
        mutated_gene=target_mutation,
        candidate_genes=candidates,
        max_papers_per_gene=max_papers_per_gene,
        email=email,
        api_key=api_key,
    )
    sections.append(literature)
    literature_scores = _parse_literature_scores(literature)

    ppi = analyze_ppi_network_for_sl(
        target_gene=target_mutation,
        candidate_genes=candidates,
        required_score=required_score,
    )
    sections.append(ppi)
    ppi_verdicts = _parse_ppi_classification(ppi)

    # Falsification stage: run it here rather than trusting the caller to remember, so a confounded
    # candidate can never reach a Go recommendation.
    confounder_verdicts = {}
    if check_confounders:
        confounders = check_dependency_confounders(
            target_mutation=target_mutation,
            candidate_genes=candidates,
            cancer_type=cancer_type,
            data_lake_path=data_lake_path,
            mutation_csv_path=mutation_csv_path,
        )
        sections.append(confounders)
        if not confounders.startswith("FAILURE"):
            current = None
            for line in confounders.splitlines():
                if line.startswith("### "):
                    current = line[4:].strip()
                elif current and line.strip().startswith("VERDICT:"):
                    confounder_verdicts[current] = line.split("VERDICT:", 1)[1].strip()
                    current = None

    # --- integration ---------------------------------------------------------------------------
    dossier = [
        "=" * 78,
        f"EVIDENCE DOSSIER - {target_mutation.upper()}-driven synthetic lethality in {cancer_type}",
        f"Generated {datetime.now(tz=_UTC).strftime('%Y-%m-%d %H:%M UTC')}",
        "=" * 78,
    ]

    ranked = []
    for gene in candidates:
        stats_row = discovery_stats.get(gene, {})
        literature_row = literature_scores.get(gene, {"score": 0, "interpretation": "not evaluated", "refuted": False})
        ppi_verdict = ppi_verdicts.get(gene, "UNKNOWN")

        effect = stats_row.get("effect_difference", 0.0)
        selectivity = stats_row.get("pct_mutant_dependent", 0.0) - 0.0
        pan_essentiality = stats_row.get("pct_all_lines_dependent", 0.0)

        # Deterministic integrated score: directness-weighted, not a rank average.
        statistical_component = min(45, abs(effect) * 60)
        selectivity_component = min(20, selectivity * 0.2)
        # Novelty must not be punished like refutation: an undocumented pair starts from a neutral
        # baseline, and only explicitly refuting literature subtracts (via the penalty below).
        literature_component = min(25, 8 + literature_row["score"] * 0.2)
        ppi_component = {
            "DIRECT-EXPERIMENTAL": 10,
            "DIRECT-HIGH-CONFIDENCE": 7,
            "DIRECT-WEAK": 4,
            "PATHWAY-PROXIMAL": 6,
            "WEAKLY-PROXIMAL": 3,
            "DISTANT": 0,
        }.get(ppi_verdict, 0)
        penalty = 0
        contradictions = []
        if pan_essentiality >= 70:
            penalty += 25
            contradictions.append(
                f"depleted in {pan_essentiality:.0f}% of all screened lines - largely a common essential gene; "
                "the genotype effect is a modulation of a general dependency and the therapeutic window is doubtful"
            )
        elif pan_essentiality >= 50:
            penalty += 15
            contradictions.append(
                f"depleted in {pan_essentiality:.0f}% of all screened lines - partly a common essential gene, "
                "so the therapeutic window in normal tissue is uncertain"
            )
        if literature_row.get("refuted"):
            penalty += 15
            contradictions.append("published refuting / non-replication statements were detected")
        if stats_row.get("q_value", 1) > 0.1:
            penalty += 5
            contradictions.append(
                f"FDR q={stats_row.get('q_value', float('nan')):.3f} - only weak multiple-testing support"
            )
        if ppi_verdict == "DISTANT":
            contradictions.append("no protein-network proximity to the driver; mechanism unexplained")

        confounder_verdict = confounder_verdicts.get(gene, "")
        confounded = confounder_verdict.startswith("CONFOUNDED")
        if confounded:
            penalty += 30
            contradictions.append(f"confounder analysis: {confounder_verdict}")

        total = max(0, statistical_component + selectivity_component + literature_component + ppi_component - penalty)

        # A confounded dependency is not about the driver at all, so it cannot earn a Go regardless
        # of how strong the statistics look.
        if confounded:
            grade, recommendation = "D", "NO-GO (confounded - the dependency is explained by something else)"
        elif total >= 70 and not literature_row.get("refuted"):
            grade, recommendation = "A", "GO"
        elif total >= 50:
            grade, recommendation = "B", "GO (with the confirmatory experiment below)"
        elif total >= 30:
            grade, recommendation = "C", "HOLD"
        else:
            grade, recommendation = "D", "NO-GO"

        ranked.append(
            {
                "gene": gene,
                "total": total,
                "grade": grade,
                "recommendation": recommendation,
                "stats": stats_row,
                "literature": literature_row,
                "ppi": ppi_verdict,
                "confounder": confounder_verdict,
                "contradictions": contradictions,
                "components": {
                    "statistical": round(statistical_component, 1),
                    "selectivity": round(selectivity_component, 1),
                    "literature": round(literature_component, 1),
                    "ppi": ppi_component,
                    "penalty": -penalty,
                },
            }
        )

    ranked.sort(key=lambda r: r["total"], reverse=True)

    dossier.append("")
    dossier.append(f"{'rank':<6}{'gene':<12}{'score':>7}{'grade':>7}  recommendation")
    dossier.append("-" * 78)
    for index, item in enumerate(ranked, start=1):
        dossier.append(f"{index:<6}{item['gene']:<12}{item['total']:>7.1f}{item['grade']:>7}  {item['recommendation']}")

    for index, item in enumerate(ranked, start=1):
        stats_row = item["stats"]
        gene = item["gene"]
        dossier.append("")
        dossier.append("-" * 78)
        dossier.append(f"CANDIDATE {index}: {target_mutation.upper()} - {gene}  ({cancer_type})")
        dossier.append("-" * 78)
        dossier.append(
            f"  Conclusion: {gene} is a candidate synthetic-lethal partner of mutant {target_mutation.upper()} "
            f"in {cancer_type} with confidence grade {item['grade']} "
            f"(integrated score {item['total']:.1f}/100)."
        )
        dossier.append(f"  Score composition: {json.dumps(item['components'])}")
        dossier.append("")
        dossier.append("  SUPPORTING EVIDENCE")
        if stats_row:
            dossier.append(
                f"    [DepMap CRISPR] gene effect {stats_row['mutant_mean_effect']:.3f} (mutant) vs "
                f"{stats_row['wildtype_mean_effect']:.3f} (wild type); delta={stats_row['effect_difference']:.3f}, "
                f"Cohen's d={stats_row['cohens_d']:.2f}, p={stats_row['p_value']:.2e}, q={stats_row['q_value']:.3f}; "
                f"{stats_row['pct_mutant_dependent']:.0f}% of mutant lines are dependent."
            )
        dossier.append(
            f"    [PubMed] support score {item['literature']['score']}/100 - {item['literature']['interpretation']}"
        )
        dossier.append(f"    [STRING PPI] {item['ppi']}")
        if item["confounder"]:
            dossier.append(f"    [Confounder check] {item['confounder']}")
        dossier.append("")
        dossier.append("  CONTRADICTING EVIDENCE / FAILURE MODES")
        if item["contradictions"]:
            for contradiction in item["contradictions"]:
                dossier.append(f"    - {contradiction}")
        else:
            dossier.append("    - No explicit contradiction detected; the main risk remains the small wild-type group.")
        dossier.append("")
        dossier.append("  MINIMAL VALIDATION EXPERIMENT")
        dossier.append(
            f"    1. Isogenic pair: knock out {gene} (2 independent sgRNAs) in a {target_mutation.upper()}-mutant "
            "line and its wild-type counterpart; read out 10-day viability plus caspase-3/7."
        )
        dossier.append(
            f"    2. Rescue: re-express sgRNA-resistant {gene} cDNA; loss of the differential effect confirms "
            "on-target action."
        )
        dossier.append(
            f"    3. Genotype panel: repeat in >=3 mutant and >=3 wild-type lines and test the "
            f"genotype x {gene}-knockout interaction term (two-way ANOVA)."
        )
        dossier.append(
            "    4. Therapeutic window: repeat in a non-transformed line (e.g. HPNE/HPDE) - a comparable "
            "dependency there is a no-go signal."
        )
        dossier.append("    Decision criterion: >=2-fold viability difference (mutant vs wild type), p<0.05, rescued.")
        dossier.append("")
        dossier.append(f"  RECOMMENDATION: {item['recommendation']}")

    dossier.append("")
    dossier.append("UNRESOLVED QUESTIONS / DATA GAPS")
    dossier.append(
        "  - All statistical evidence derives from a single DepMap release; results were not replicated in an "
        "independent screen (e.g. Sanger Project Score) and must not be counted as independent confirmation."
    )
    dossier.append(
        "  - Single-gene knockout dependency is a proxy for synthetic lethality; no combinatorial CRISPR data is used."
    )
    dossier.append(
        "  - Literature scores measure documentation, not truth: an unexplored pair scores low yet may still be real."
    )
    dossier.append(
        "  - Cell-line dependency does not establish a therapeutic window in patients; normal-tissue toxicity is untested."
    )

    # The dossier goes first: agent frameworks crop tool output to the first ~10,000 characters,
    # and the Go/Hold/No-go recommendation is precisely the part that must survive that crop.
    output = ["\n".join(dossier)]
    if include_stage_logs:
        output.append("\n\n".join(["=" * 78, "PER-STAGE RESEARCH LOGS (raw tool output)", *sections]))
    else:
        output.append(
            "Per-stage research logs were omitted to keep this output readable. Re-run with "
            "include_stage_logs=True, or call discover_synthetic_lethal_candidates / "
            "validate_sl_candidates_with_pubmed / analyze_ppi_network_for_sl / "
            "check_dependency_confounders individually, to see the raw evidence."
        )
    return "\n\n".join(output)


# ---------------------------------------------------------------------------
# Tool 6: ovarian cancer mutation-stratified dependency analysis
# ---------------------------------------------------------------------------
def stratify_ovarian_cancer_dependency_by_mutation(
    target_mutation: str,
    histology: str | None = None,
    stratify_by: str = "mutation",
    data_lake_path: str | None = None,
    mutation_csv_path: str | None = None,
    candidate_genes=None,
    restrict_to_synlethdb_partners: bool = False,
    top_n: int = 25,
    p_threshold: float = 0.05,
    fdr_threshold: float = 0.25,
    min_effect_difference: float = -0.2,
    max_mutant_mean_effect: float = -0.3,
    exclude_pan_essential: bool = True,
    output_csv_path: str | None = None,
    plot_output_prefix: str | None = None,
) -> str:
    """Find genes whose CRISPR knockout hits altered ovarian cancer lines harder than unaltered lines.

    DepMap ovarian cell lines (OncotreeLineage "Ovary/Fallopian Tube", optionally restricted to one
    histology such as high-grade serous) are split into an altered group and a control group by the
    genotype of ``target_mutation``: somatic mutation (default), copy-number amplification, or deep
    deletion. Amplification matters in ovarian cancer specifically because CCNE1, one of the central
    HGSOC drivers, is an amplification event that mutation calls cannot see at all. Every tested gene
    is compared between the two groups with a Welch t-test
    (unequal variances, which is the right test here because the mutant group is usually much
    smaller than the control group) on the Chronos gene-effect score, where more negative means
    stronger dependency. The reported p-value is the one-sided p for the directional hypothesis
    "altered lines are MORE depleted", multiplicity-corrected with Benjamini-Hochberg across the
    genes actually tested. Surviving genes are annotated with SynLethDB evidence so that a hit
    already recorded as a synthetic-lethal partner of the driver can be told apart from a novel one.

    Parameters
    ----------
    target_mutation : str
        HUGO symbol of the gene whose alteration defines the two groups, e.g. "BRCA1", "ARID1A",
        "TP53", "CCNE1".
    histology : str, optional
        Restrict the cohort to one ovarian histology. Accepts "HGSOC" / "high-grade serous",
        "clear cell", "endometrioid", "mucinous", "serous", or any substring of OncotreeSubtype.
        Default: all ovarian tumour lines.
    stratify_by : str, optional
        Genotype event that splits the cohort: "mutation" (default), "amplification" (use for
        CCNE1, MYC, ERBB2 and other amplification drivers) or "deletion" (homozygous loss, e.g.
        PTEN, RB1). Copy-number calls come from a local DepMap OmicsCNGene.csv when present,
        otherwise from the cBioPortal CCLE discrete CNA profile - no extra download is required.
    data_lake_path : str, optional
        Directory holding DepMap_CRISPRGeneEffect.csv, DepMap_Model.csv and (optionally)
        synlethdb_human_sl.parquet, DepMap_OmicsCNGene.csv (default: "./data/biomni_data/data_lake").
    mutation_csv_path : str, optional
        Genotype table overriding the default call source, for any ``stratify_by`` mode. Columns:
        ModelID, HugoSymbol[, ProteinChange] for mutations, or ModelID, HugoSymbol[, Alteration]
        for copy number, where Alteration contains "AMP" or "DEL". Without it, a local DepMap file
        is used when present, otherwise calls come from the cBioPortal CCLE study.
    candidate_genes : list[str] | str, optional
        Restrict testing to these genes (list, or comma/whitespace separated string). A focused
        gene set makes the FDR correction far less punishing than a genome-wide scan.
    restrict_to_synlethdb_partners : bool, optional
        Test only genes recorded in SynLethDB as synthetic-lethal partners of ``target_mutation``
        (default: False). Combines with ``candidate_genes`` as an intersection.
    top_n : int, optional
        Number of top candidates to print in detail (default: 25).
    p_threshold : float, optional
        One-sided Welch t-test p-value cutoff (default: 0.05).
    fdr_threshold : float, optional
        Benjamini-Hochberg q-value cutoff (default: 0.25). Set to 1.0 to disable FDR filtering.
    min_effect_difference : float, optional
        Required (altered-group mean - control-group mean) gene effect; must be negative (default: -0.2).
    max_mutant_mean_effect : float, optional
        The mutant group mean gene effect must be below this value, so that statistically
        significant but biologically trivial differences are dropped (default: -0.3).
    exclude_pan_essential : bool, optional
        Drop common-essential genes depleted in >80% of all screened lines (default: True).
    output_csv_path : str, optional
        If given, the full ranked table (all tested genes) is written to this CSV path.
    plot_output_prefix : str, optional
        If given, writes "<prefix>_volcano.png" (all tested genes) and "<prefix>_boxplot.png"
        (per-cell-line dependency of the top candidates in both groups) and reports their paths.

    Returns
    -------
    str
        A research log with the ovarian cohort composition, altered/control group membership,
        the ranked differential-dependency table with SynLethDB annotation, QC warnings and
        full provenance.

    """
    import numpy as np
    import pandas as pd
    from scipy import stats

    driver = str(target_mutation).strip().upper()
    mode = _STRATIFY_ALIASES.get(str(stratify_by).strip().lower())
    if mode is None:
        return (
            f"FAILURE: stratify_by='{stratify_by}' is not recognised. Use one of "
            f"{sorted(set(_STRATIFY_ALIASES.values()))}."
        )
    altered_label, control_label, event_name, (short_alt, short_ctrl) = _STRATIFY_MODES[mode]

    log = [
        "=" * 78,
        f"OVARIAN CANCER {event_name.upper()}-STRATIFIED DEPENDENCY - {driver}"
        + (f" ({histology})" if histology else " (all ovarian histologies)"),
        f"Run at {datetime.now(tz=_UTC).strftime('%Y-%m-%d %H:%M UTC')}",
        "=" * 78,
        "",
    ]

    try:
        groups = _ovarian_groups(driver, histology, mode, data_lake_path, mutation_csv_path)
    except (FileNotFoundError, ValueError, RuntimeError) as e:
        return f"FAILURE: {e}"

    bundle = groups["bundle"]
    gene_effect = bundle["gene_effect"]
    cohort = groups["cohort"]
    cohort_note = groups["cohort_note"]
    genotype_source = groups["genotype_source"]
    mutant_ids = groups["mutant_ids"]
    wildtype_ids = groups["wildtype_ids"]
    unknown_ids = groups["unknown_ids"]

    # --- Step 1: ovarian cohort -------------------------------------------------------------
    subtype_counts = cohort["OncotreeSubtype"].value_counts().to_dict()
    log.append("STEP 1 | Ovarian cohort definition")
    log.append(f"  {cohort_note}")
    log.append(f"  Lines with CRISPR gene-effect data: {len(cohort)}")
    log.append(f"  Histology composition: {subtype_counts}")

    # --- Step 2: genotype stratification ------------------------------------------------------
    log.append("")
    log.append(f"STEP 2 | Genotype stratification ({event_name})")
    log.append(f"  {event_name.capitalize()} source: {genotype_source}")
    log.append(f"  {driver}-{altered_label} lines: {len(mutant_ids)}")
    log.append(f"  {driver}-{control_label} lines: {len(wildtype_ids)}")
    log.append(f"  Not profiled (excluded) : {len(unknown_ids)}")
    variant_counts = cohort.loc[cohort["MutationStatus"] == "MUT", "Variant"].value_counts().head(6).to_dict()
    if variant_counts:
        log.append(f"  Alteration spectrum (top 6): {variant_counts}")
    log.append(
        f"  {altered_label}: "
        f"{', '.join(cohort.loc[cohort['MutationStatus'] == 'MUT', 'StrippedCellLineName'].tolist())}"
    )
    log.append(
        f"  {control_label}: "
        f"{', '.join(cohort.loc[cohort['MutationStatus'] == 'WT', 'StrippedCellLineName'].tolist())}"
    )

    if len(mutant_ids) < 3 or len(wildtype_ids) < 3:
        log.append("")
        log.append(
            f"FAILURE: insufficient group sizes ({altered_label.lower()} n={len(mutant_ids)}, "
            f"{control_label.lower()} n={len(wildtype_ids)}); "
            "at least 3 per group are required for a Welch t-test. Drop the `histology` filter, or use "
            "discover_synthetic_lethal_candidates with cancer_type='pan-cancer' for a larger cohort."
        )
        return "\n".join(log)

    # --- Step 3: gene set to test -------------------------------------------------------------
    synlethdb = _load_synlethdb(bundle["data_lake_path"])
    partners = _synlethdb_partners(driver, bundle["data_lake_path"])

    tested_columns = list(gene_effect.columns)
    restriction_notes = []
    requested = _parse_gene_list(candidate_genes)
    if requested:
        keep = {bundle["gene_columns"][g] for g in requested if g in bundle["gene_columns"]}
        missing = [g for g in requested if g not in bundle["gene_columns"]]
        tested_columns = [c for c in tested_columns if c in keep]
        restriction_notes.append(f"candidate_genes: {len(requested)} requested, {len(keep)} in the CRISPR library")
        if missing:
            restriction_notes.append(f"  not screened / unknown symbol: {', '.join(missing[:20])}")
    if restrict_to_synlethdb_partners:
        if synlethdb is None:
            return (
                "FAILURE: restrict_to_synlethdb_partners=True but synlethdb_human_sl.parquet was not found in "
                f"{bundle['data_lake_path']}. Provide the file or set the flag to False."
            )
        keep = {bundle["gene_columns"][g] for g in partners if g in bundle["gene_columns"]}
        tested_columns = [c for c in tested_columns if c in keep]
        restriction_notes.append(
            f"SynLethDB partners of {driver}: {len(partners)} recorded, {len(keep)} screened in DepMap"
        )
    if requested or restrict_to_synlethdb_partners:
        # Keep the driver itself in the tested set, otherwise the self-dependency control below
        # has nothing to check the mutation calls against.
        driver_column = bundle["gene_columns"].get(driver)
        if driver_column and driver_column not in tested_columns:
            tested_columns.append(driver_column)
    if len(tested_columns) == 0:
        return (
            f"FAILURE: the gene restriction left no testable gene. {'; '.join(restriction_notes)}"
        )

    # --- Step 4: differential dependency testing ----------------------------------------------
    mutant_matrix = gene_effect.loc[mutant_ids, tested_columns]
    wildtype_matrix = gene_effect.loc[wildtype_ids, tested_columns]
    usable = mutant_matrix.columns[(mutant_matrix.notna().sum() >= 3) & (wildtype_matrix.notna().sum() >= 3)]
    mutant_matrix = mutant_matrix[usable]
    wildtype_matrix = wildtype_matrix[usable]

    tstat, p_two_sided = stats.ttest_ind(
        mutant_matrix.values, wildtype_matrix.values, axis=0, equal_var=False, nan_policy="omit"
    )
    tstat = np.asarray(tstat, dtype=float)
    p_two_sided = np.asarray(p_two_sided, dtype=float)
    p_two_sided = np.where(np.isfinite(p_two_sided), p_two_sided, 1.0)
    # Directional hypothesis: the mutant group is MORE depleted, i.e. t < 0.
    p_one_sided = np.where(tstat < 0, p_two_sided / 2.0, 1.0 - p_two_sided / 2.0)

    mutant_mean = np.asarray(mutant_matrix.mean())
    wildtype_mean = np.asarray(wildtype_matrix.mean())
    pooled_sd = np.sqrt((np.asarray(mutant_matrix.std(ddof=1)) ** 2 + np.asarray(wildtype_matrix.std(ddof=1)) ** 2) / 2)
    with np.errstate(divide="ignore", invalid="ignore"):
        cohens_d = np.where(pooled_sd > 0, (mutant_mean - wildtype_mean) / pooled_sd, np.nan)

    # Pan-essentiality is judged across the whole screen, not just the ovarian cohort.
    pan_fraction = (gene_effect[usable] < DEPLETION_THRESHOLD).sum() / gene_effect[usable].notna().sum()
    mutant_dependent = (mutant_matrix < DEPLETION_THRESHOLD).sum() / mutant_matrix.notna().sum()
    wildtype_dependent = (wildtype_matrix < DEPLETION_THRESHOLD).sum() / wildtype_matrix.notna().sum()

    results = pd.DataFrame(
        {
            "gene": [c.split(" (")[0] for c in usable],
            "n_mutant": np.asarray(mutant_matrix.notna().sum()),
            "n_wildtype": np.asarray(wildtype_matrix.notna().sum()),
            "mutant_mean_effect": mutant_mean,
            "wildtype_mean_effect": wildtype_mean,
            "effect_difference": mutant_mean - wildtype_mean,
            "cohens_d": cohens_d,
            "t_statistic": tstat,
            "p_one_sided": p_one_sided,
            "p_two_sided": p_two_sided,
            "q_value": _benjamini_hochberg(p_one_sided),
            "pct_mutant_dependent": np.asarray(mutant_dependent) * 100,
            "pct_wildtype_dependent": np.asarray(wildtype_dependent) * 100,
            "pct_all_lines_dependent": np.asarray(pan_fraction) * 100,
        }
    )
    results["pan_essential"] = results["pct_all_lines_dependent"] >= PAN_ESSENTIAL_FRACTION * 100
    results["synlethdb_known"] = results["gene"].str.upper().isin(partners)
    results["synlethdb_evidence"] = results["gene"].str.upper().map(
        lambda g: ";".join(partners[g]["sources"]) if g in partners else ""
    )
    results["synlethdb_tier"] = results["gene"].str.upper().map(
        lambda g: _SYNLETHDB_TIER_NAME[partners[g]["tier"]] if g in partners else ""
    )
    # Rank over all tested genes, used below to state where known SL partners land.
    results["rank_by_delta"] = results["effect_difference"].rank(method="min").astype(int)

    log.append("")
    log.append(
        f"STEP 3 | Differential dependency testing (Welch t-test, {altered_label.lower()} vs {control_label.lower()})"
    )
    for note in restriction_notes:
        log.append(f"  {note}")
    log.append(f"  Genes tested: {len(results)} (of {gene_effect.shape[1]} in the CRISPR library)")
    log.append(
        f"  One-sided p < {p_threshold} ({altered_label.lower()} more depleted): "
        f"{(results['p_one_sided'] < p_threshold).sum()}"
    )
    log.append(f"  BH q < {fdr_threshold}: {(results['q_value'] < fdr_threshold).sum()}")

    # --- Step 5: filtering / QC ---------------------------------------------------------------
    selected = results[
        (results["p_one_sided"] < p_threshold)
        & (results["q_value"] < fdr_threshold)
        & (results["effect_difference"] <= min_effect_difference)
        & (results["mutant_mean_effect"] <= max_mutant_mean_effect)
    ].copy()
    n_before_pan = len(selected)
    if exclude_pan_essential:
        selected = selected[~selected["pan_essential"]]
    selected = selected.sort_values(["effect_difference", "p_one_sided"]).reset_index(drop=True)
    selected.insert(0, "rank", np.arange(1, len(selected) + 1))

    log.append("")
    log.append("STEP 4 | Filtering and quality control")
    log.append(f"  After significance + effect-size filters: {n_before_pan}")
    log.append(
        f"  Pan-essential removed (depleted in >{PAN_ESSENTIAL_FRACTION:.0%} of all "
        f"{gene_effect.shape[0]} screened lines): {n_before_pan - len(selected)}"
    )
    log.append(f"  Final candidate count: {len(selected)}")

    control = results[results["gene"].str.upper() == driver]
    if len(control):
        row = control.iloc[0]
        right_direction = row["effect_difference"] <= min_effect_difference
        if right_direction and row["p_one_sided"] < p_threshold:
            verdict = "PASS"
        elif right_direction:
            verdict = "WEAK"
        else:
            verdict = "FAIL"
        log.append(
            f"  Self-dependency control ({verdict}): {driver} itself delta={row['effect_difference']:.3f}, "
            f"p1={row['p_one_sided']:.2e} ({altered_label} {row['mutant_mean_effect']:.3f} vs "
            f"{control_label} {row['wildtype_mean_effect']:.3f})."
        )
        if verdict == "WEAK":
            log.append(
                "    The self-dependency has the right direction and size but misses the p cutoff - that is a "
                "power limit of these group sizes, not a contradiction of the genotype calls."
            )
        elif verdict == "FAIL" and mode == "amplification":
            log.append(
                "    WARNING: an amplified oncogene should be a self-dependency in its own amplified lines. "
                "A FAIL here means the copy-number calls or the cohort are mis-specified, and every candidate "
                "below should be treated as unreliable until that is resolved."
            )
        elif verdict == "FAIL":
            log.append(
                "    Note: a FAIL is expected for a tumour suppressor (loss-of-function drivers are not "
                "self-dependencies) but is a red flag for an oncogene / amplification driver."
            )
    else:
        log.append(f"  Self-dependency control unavailable: {driver} was not among the tested genes.")

    if synlethdb is None:
        log.append("  SynLethDB annotation unavailable: synlethdb_human_sl.parquet not found in the data lake.")
    else:
        tested_partners = results[results["synlethdb_known"]]
        selected_partners = selected[selected["synlethdb_known"]] if len(selected) else selected
        log.append(
            f"  SynLethDB: {len(partners)} known SL partners of {driver}, {len(tested_partners)} of them tested here, "
            f"{len(selected_partners)} in the final candidate list."
        )
        if len(tested_partners):
            best = tested_partners.nsmallest(5, "effect_difference")
            log.append(
                "    Best-ranking known partners (rank of "
                f"{len(results)} tested, by delta): "
                + ", ".join(
                    f"{r['gene']}#{r['rank_by_delta']} (delta={r['effect_difference']:.3f}, p1={r['p_one_sided']:.1e})"
                    for _, r in best.iterrows()
                )
            )
            log.append(
                "    Known partners are a calibration read-out, not a filter: if none of them rank anywhere "
                "near the top, the cohort is underpowered and the novel hits below are unreliable."
            )

    warnings = []
    if len(mutant_ids) < 8 or len(wildtype_ids) < 8:
        warnings.append(
            f"Low power: {altered_label.lower()} n={len(mutant_ids)}, {control_label.lower()} "
            f"n={len(wildtype_ids)}. A single outlier line can "
            "create or destroy a candidate; treat p-values as descriptive and replicate in a larger cohort."
        )
    mutant_fraction = len(mutant_ids) / (len(mutant_ids) + len(wildtype_ids))
    if mutant_fraction > 0.9 or mutant_fraction < 0.1:
        warnings.append(
            f"Unbalanced groups ({mutant_fraction:.0%} {altered_label.lower()}). For near-universal ovarian "
            "drivers such as TP53 in "
            "high-grade serous disease the wild-type group is not a comparable control - the few TP53-wild-type "
            "ovarian lines are usually a different histology altogether."
        )
    if len(unknown_ids):
        warnings.append(
            f"{len(unknown_ids)} lines had no {event_name} call and were excluded rather than assumed "
            f"{control_label.lower()}."
        )
    if len(subtype_counts) > 1:
        warnings.append(
            f"The cohort mixes {len(subtype_counts)} histologies {subtype_counts}; histology, not the mutation, "
            "may drive some contrasts. Re-run with `histology` set to test within one subtype."
        )
    if synlethdb is not None and synlethdb["constant_score"]:
        warnings.append(
            "SynLethDB `score` is constant (1.0) in this snapshot, so it cannot rank confidence; the evidence "
            "tier printed per candidate comes from the `source` field instead, and 'Computational Prediction' or "
            "'Text Mining' support is not experimental validation."
        )
    if mode == "mutation":
        warnings.append(
            "Mutation calls are presence/absence of any protein-altering variant - they do not distinguish "
            "loss-of-function from VUS, monoallelic from biallelic loss, or mutation from BRCA1 promoter "
            "methylation, which is a common HR-deficiency mechanism in ovarian cancer and leaves such lines in "
            "the wild-type group."
        )
    else:
        warnings.append(
            f"{event_name.capitalize()} is a DNA-level call: it does not prove the gene is over- or "
            "under-expressed in these lines. Confirm with the DepMap expression matrix "
            "(check_dependency_confounders does this) before treating the copy-number group as a functional one."
        )
        if "cBioPortal" in genotype_source:
            warnings.append(
                "Copy-number calls come from the cBioPortal CCLE 2019 snapshot as discrete GISTIC-style values "
                "(-2/0/2 only, matched by cell line name), not from the current DepMap release. Shallow "
                "single-copy events are invisible, and a few lines may fail name matching. Put a current "
                "OmicsCNGene.csv in the data lake as DepMap_OmicsCNGene.csv to use continuous DepMap calls instead."
            )
    warnings.append(
        "Single-gene knockout dependency in cell lines is correlative evidence for a genotype-selective "
        "vulnerability, not a validated synthetic lethal interaction."
    )

    log.append("")
    log.append("QC WARNINGS (falsification checklist)")
    for warning in warnings:
        log.append(f"  - {warning}")

    # --- Step 6: report -----------------------------------------------------------------------
    log.append("")
    log.append(f"TOP {min(top_n, len(selected))} {altered_label}-SELECTIVE DEPENDENCIES")
    log.append(
        f"{'rank':<5}{'gene':<12}{short_alt + ' mean':>10}{short_ctrl + ' mean':>10}"
        f"{'delta':>9}{'d':>7}"
        f"{'p(1-sided)':>12}{'q':>9}{'%alt dep':>10}{'%all dep':>10}  SynLethDB"
    )
    log.append("-" * 105)
    for _, row in selected.head(top_n).iterrows():
        sl_note = f"known: {row['synlethdb_tier']}" if row["synlethdb_known"] else "-"
        log.append(
            f"{int(row['rank']):<5}{row['gene']:<12}{row['mutant_mean_effect']:>10.3f}"
            f"{row['wildtype_mean_effect']:>10.3f}{row['effect_difference']:>9.3f}{row['cohens_d']:>7.2f}"
            f"{row['p_one_sided']:>12.2e}{row['q_value']:>9.3f}{row['pct_mutant_dependent']:>9.0f}%"
            f"{row['pct_all_lines_dependent']:>9.0f}%  {sl_note}"
        )

    if len(selected) == 0:
        # An empty table is the normal outcome of a genome-wide FDR gate on a cohort this small, so
        # report what the run actually saw instead of nothing. These are NOT candidates.
        runners = results[
            (results["p_one_sided"] < p_threshold)
            & (results["effect_difference"] <= min_effect_difference)
            & (results["mutant_mean_effect"] <= max_mutant_mean_effect)
        ]
        if exclude_pan_essential:
            runners = runners[~runners["pan_essential"]]
        runners = runners.sort_values(["effect_difference", "p_one_sided"]).head(top_n)
        log.append("")
        log.append(
            f"NO GENE SURVIVED THE FDR GATE (q < {fdr_threshold}). Nominally significant runners-up "
            f"(p < {p_threshold} but NOT significant after correction for {len(results)} tests):"
        )
        for _, row in runners.iterrows():
            sl_note = f"  [SynLethDB known: {row['synlethdb_tier']}]" if row["synlethdb_known"] else ""
            log.append(
                f"  {row['gene']:<12}{short_alt} {row['mutant_mean_effect']:>7.3f}  "
                f"{short_ctrl} {row['wildtype_mean_effect']:>7.3f}  "
                f"delta {row['effect_difference']:>7.3f}  p1 {row['p_one_sided']:.2e}  q {row['q_value']:.3f}{sl_note}"
            )
        advice = (
            "widen the cohort (drop `histology`, or use discover_synthetic_lethal_candidates with "
            "cancer_type='pan-cancer') - the gene set is already as small as this analysis can make it"
            if restrict_to_synlethdb_partners
            else "pass `candidate_genes` (a pathway or druggable gene set) or `restrict_to_synlethdb_partners=True`"
        )
        log.append(
            "  These are hypotheses at nominal significance only. To get a calibrated candidate list, shrink the "
            f"multiple-testing burden rather than raising the FDR cutoff: {advice}. An altered group of this size "
            "cannot clear BH correction over a large gene set, which is a power limit of the cohort, not evidence "
            "of no effect."
        )

    candidate_list = [g for g in selected["gene"].head(top_n).tolist() if g.upper() != driver]
    log.append("")
    log.append(f"CANDIDATE_GENES: {', '.join(candidate_list)}")
    log.append(
        "  (pass this list to validate_sl_candidates_with_pubmed, analyze_ppi_network_for_sl and "
        "check_dependency_confounders for literature, network and confounder evidence)"
    )

    if output_csv_path:
        os.makedirs(os.path.dirname(os.path.abspath(output_csv_path)), exist_ok=True)
        results.sort_values(["effect_difference", "p_one_sided"]).to_csv(output_csv_path, index=False)
        log.append(f"  Full table of all {len(results)} tested genes written to {output_csv_path}")

    if plot_output_prefix:
        log.append("")
        log.append("FIGURES")
        cohort_title = (
            "ovarian" + (f" ({histology})" if histology else "")
            + f": {altered_label} n={len(mutant_ids)} vs {control_label} n={len(wildtype_ids)}"
        )
        volcano_path = f"{plot_output_prefix}_volcano.png"
        try:
            counts = _draw_volcano(
                results, volcano_path, f"{driver} {event_name} - {cohort_title}", p_threshold, fdr_threshold,
                min_effect_difference, exclude_pan_essential, 15,
            )
            log.append(
                f"  Volcano plot: {volcano_path} ({counts['n_total']} genes, {counts['n_selected']} significant, "
                f"{counts['n_nominal']} nominal only)"
            )
        except Exception as e:  # a plotting failure must never lose the analysis above
            log.append(f"  Volcano plot FAILED: {e}")

        if len(selected):
            plot_genes = selected["gene"].head(8).tolist()
            plot_note = "top candidates"
        else:
            fallback = results[
                (results["p_one_sided"] < p_threshold)
                & (results["effect_difference"] <= min_effect_difference)
                & (results["mutant_mean_effect"] <= max_mutant_mean_effect)
            ]
            if exclude_pan_essential:
                fallback = fallback[~fallback["pan_essential"]]
            plot_genes = fallback.nsmallest(8, "effect_difference")["gene"].tolist()
            plot_note = "nominally significant runners-up - no gene cleared the FDR gate"
        boxplot_path = f"{plot_output_prefix}_boxplot.png"
        if plot_genes:
            try:
                _draw_dependency_boxplot(
                    bundle, cohort, mutant_ids, wildtype_ids, plot_genes,
                    (altered_label, control_label), boxplot_path,
                    f"CRISPR dependency by {driver} {event_name} status - {cohort_title}",
                )
                log.append(f"  Boxplot: {boxplot_path} ({plot_note}: {', '.join(plot_genes)})")
            except Exception as e:
                log.append(f"  Boxplot FAILED: {e}")
        else:
            log.append("  Boxplot skipped: no gene passed even the nominal significance and effect-size filters.")

    log.append("")
    log.append("PROVENANCE")
    for item in bundle["provenance"]:
        log.append(f"  - {item}")
    log.append(f"  - {event_name.capitalize()} calls: {genotype_source}")
    if synlethdb is not None:
        log.append(f"  - {synlethdb['provenance']}")
    log.append(
        f"  - Statistics: Welch two-sample t-test (equal_var=False) on Chronos gene effect, one-sided for "
        f"'{altered_label.lower()} more depleted', BH-FDR across {len(results)} tested genes; depletion threshold "
        f"{DEPLETION_THRESHOLD}, pan-essential cutoff {PAN_ESSENTIAL_FRACTION:.0%}"
    )
    return "\n".join(log)


# ---------------------------------------------------------------------------
# Visualisation
# ---------------------------------------------------------------------------
# Group colours are kept fixed across every figure so that "altered" always reads as the same
# colour whether the stratifying event is a mutation, an amplification or a deletion.
ALTERED_COLOR = "#d1495b"
CONTROL_COLOR = "#00798c"


def _apply_plot_style():
    """Headless-safe matplotlib/seaborn setup. Must run before pyplot is used."""
    import matplotlib

    matplotlib.use("Agg")
    import seaborn as sns

    sns.set_theme(style="whitegrid", context="notebook")
    return sns


def _ovarian_groups(
    driver: str,
    histology: str | None,
    mode: str,
    data_lake_path: str | None,
    genotype_csv_path: str | None,
) -> dict:
    """Build the ovarian cohort and split it into altered / control groups.

    Shared by the analysis tool and the plotting tools so that a figure always shows exactly the
    cell lines the statistics were computed on. Raises ValueError / RuntimeError with a
    human-readable message when the cohort cannot be built.
    """
    bundle = _load_depmap(data_lake_path)
    gene_effect = bundle["gene_effect"]

    cohort, cohort_note = _select_ovarian_models(bundle["model"], histology)
    cohort = cohort[cohort["ModelID"].isin(gene_effect.index)]
    if len(cohort) == 0:
        subtypes = sorted(_select_ovarian_models(bundle["model"])[0]["OncotreeSubtype"].dropna().unique().tolist())
        raise ValueError(
            f"no ovarian DepMap line with CRISPR data matched histology='{histology}'. "
            f"Available ovarian subtypes: {subtypes}"
        )

    if mode == "mutation":
        cohort, genotype_source = _annotate_mutation_status(
            cohort, driver, bundle["data_lake_path"], genotype_csv_path
        )
    else:
        cohort, genotype_source = _annotate_copy_number_status(
            cohort, driver, bundle["data_lake_path"], mode, genotype_csv_path
        )

    return {
        "bundle": bundle,
        "cohort": cohort,
        "cohort_note": cohort_note,
        "genotype_source": genotype_source,
        "mutant_ids": cohort.loc[cohort["MutationStatus"] == "MUT", "ModelID"].tolist(),
        "wildtype_ids": cohort.loc[cohort["MutationStatus"] == "WT", "ModelID"].tolist(),
        "unknown_ids": cohort.loc[cohort["MutationStatus"] == "UNKNOWN", "ModelID"].tolist(),
    }


def _draw_dependency_boxplot(
    bundle: dict,
    cohort,
    mutant_ids: list,
    wildtype_ids: list,
    genes: list,
    labels: tuple,
    output_path: str,
    title: str,
) -> tuple:
    """Draw one box + strip panel per gene comparing gene effect between the two groups."""
    import numpy as np
    import pandas as pd
    from scipy import stats

    sns = _apply_plot_style()
    import matplotlib.pyplot as plt

    altered_label, control_label = labels
    gene_effect = bundle["gene_effect"]
    name_by_model = dict(zip(cohort["ModelID"], cohort["StrippedCellLineName"], strict=False))

    plotted, skipped = [], []
    records = []
    for gene in genes:
        column = bundle["gene_columns"].get(gene.upper())
        if column is None:
            skipped.append(gene)
            continue
        plotted.append(gene)
        for group, model_ids in ((altered_label, mutant_ids), (control_label, wildtype_ids)):
            values = gene_effect.loc[model_ids, column].dropna()
            for model_id, value in values.items():
                records.append(
                    {
                        "gene": gene.upper(),
                        "group": group,
                        "gene_effect": float(value),
                        "cell_line": name_by_model.get(model_id, model_id),
                    }
                )
    if not plotted:
        raise ValueError(f"none of the requested genes are in the CRISPR library: {', '.join(genes)}")

    frame = pd.DataFrame(records)
    n_columns = min(4, len(plotted))
    n_rows = int(np.ceil(len(plotted) / n_columns))
    # Shared y: the panels are all in the same gene-effect units, so a per-panel axis would make
    # a weak difference look as dramatic as a strong one.
    figure, axes = plt.subplots(
        n_rows, n_columns, figsize=(3.4 * n_columns, 3.9 * n_rows), squeeze=False, sharey=True
    )

    palette = {altered_label: ALTERED_COLOR, control_label: CONTROL_COLOR}
    for index, gene in enumerate(plotted):
        axis = axes[index // n_columns][index % n_columns]
        subset = frame[frame["gene"] == gene.upper()]
        sns.boxplot(
            data=subset, x="group", y="gene_effect", hue="group", order=[altered_label, control_label],
            palette=palette, width=0.55, showfliers=False, legend=False, ax=axis,
        )
        sns.stripplot(
            data=subset, x="group", y="gene_effect", order=[altered_label, control_label],
            color="0.2", size=4, jitter=0.18, alpha=0.75, ax=axis,
        )

        altered_values = subset.loc[subset["group"] == altered_label, "gene_effect"].to_numpy()
        control_values = subset.loc[subset["group"] == control_label, "gene_effect"].to_numpy()
        delta = float(np.mean(altered_values) - np.mean(control_values))
        tstat, p_two = stats.ttest_ind(altered_values, control_values, equal_var=False)
        p_one = p_two / 2 if tstat < 0 else 1 - p_two / 2

        # 0 = no effect, -0.5 = the depletion call used throughout this module.
        axis.axhline(0, color="0.35", linewidth=0.9)
        axis.axhline(DEPLETION_THRESHOLD, color="0.55", linewidth=0.9, linestyle="--")
        axis.set_title(
            f"{gene.upper()}\n$\\Delta$={delta:+.3f}, one-sided p={p_one:.1e}\n"
            f"n={len(altered_values)} vs {len(control_values)}",
            fontsize=10,
        )
        axis.set_xlabel("")
        axis.set_ylabel("CRISPR gene effect (Chronos)" if index % n_columns == 0 else "")
        axis.tick_params(axis="x", labelsize=9)

    for empty in range(len(plotted), n_rows * n_columns):
        axes[empty // n_columns][empty % n_columns].axis("off")

    figure.suptitle(title, fontsize=12, y=0.995)
    figure.text(
        0.5, 0.005,
        "More negative = stronger dependency. Dashed line = depletion threshold "
        f"({DEPLETION_THRESHOLD}); solid line = no effect.",
        ha="center", fontsize=8.5, color="0.35",
    )
    figure.tight_layout(rect=(0, 0.02, 1, 0.98))
    os.makedirs(os.path.dirname(os.path.abspath(output_path)), exist_ok=True)
    figure.savefig(output_path, dpi=200, bbox_inches="tight")
    plt.close(figure)
    return plotted, skipped


def _draw_volcano(
    results,
    output_path: str,
    title: str,
    p_threshold: float,
    fdr_threshold: float,
    min_effect_difference: float,
    exclude_pan_essential: bool,
    label_top: int,
) -> dict:
    """Draw effect difference vs -log10(one-sided p) over every tested gene."""
    import numpy as np

    _apply_plot_style()
    import matplotlib.pyplot as plt

    frame = results.copy()
    if "p_one_sided" not in frame.columns and "p_value" in frame.columns:
        frame["p_one_sided"] = frame["p_value"]
    for required in ("gene", "effect_difference", "p_one_sided"):
        if required not in frame.columns:
            raise ValueError(f"the results table needs a '{required}' column; found {list(frame.columns)}")
    if "q_value" not in frame.columns:
        frame["q_value"] = np.nan
    if "pan_essential" not in frame.columns:
        frame["pan_essential"] = False
    if "synlethdb_known" not in frame.columns:
        frame["synlethdb_known"] = False

    # A p of exactly 0 would be an infinite y; clip to the smallest representable p instead.
    frame["neg_log10_p"] = -np.log10(np.clip(frame["p_one_sided"].astype(float), 1e-300, 1.0))

    passes_effect = (frame["effect_difference"] <= min_effect_difference) & (frame["p_one_sided"] < p_threshold)
    passes_fdr = frame["q_value"] < fdr_threshold
    is_pan = frame["pan_essential"].fillna(False).astype(bool)
    selected = passes_effect & passes_fdr & ~(is_pan if exclude_pan_essential else False)
    nominal = passes_effect & ~selected & ~is_pan

    figure, axis = plt.subplots(figsize=(9, 7))
    background = frame[~(selected | nominal | is_pan)]
    axis.scatter(
        background["effect_difference"], background["neg_log10_p"],
        s=8, color="0.75", alpha=0.45, linewidths=0, rasterized=True, label="not significant",
    )
    pan = frame[is_pan]
    if len(pan):
        axis.scatter(
            pan["effect_difference"], pan["neg_log10_p"],
            s=14, color="0.45", alpha=0.6, marker="x", linewidths=0.8,
            label=f"pan-essential (>{PAN_ESSENTIAL_FRACTION:.0%} of all lines)",
        )
    axis.scatter(
        frame.loc[nominal, "effect_difference"], frame.loc[nominal, "neg_log10_p"],
        s=26, color="#f0a202", alpha=0.85, linewidths=0, label=f"nominal only (p < {p_threshold})",
    )
    axis.scatter(
        frame.loc[selected, "effect_difference"], frame.loc[selected, "neg_log10_p"],
        s=48, color=ALTERED_COLOR, edgecolor="black", linewidths=0.5,
        label=f"significant (q < {fdr_threshold})",
    )
    known = frame[(selected | nominal) & frame["synlethdb_known"].fillna(False).astype(bool)]
    if len(known):
        axis.scatter(
            known["effect_difference"], known["neg_log10_p"],
            s=150, facecolor="none", edgecolor="#2e294e", linewidths=1.4, label="known SynLethDB partner",
        )

    axis.axvline(0, color="0.35", linewidth=0.9)
    axis.axvline(min_effect_difference, color="0.5", linestyle="--", linewidth=0.9)
    axis.axhline(-np.log10(p_threshold), color="0.5", linestyle="--", linewidth=0.9)

    # Label the strongest hits, skipping any label that would sit on top of one already placed -
    # in a genome-wide run the significant region is dense enough that unfiltered labels overlap
    # into an unreadable block. Candidates are drawn from a wider pool so the quota is still filled.
    labelled = frame[selected | nominal]
    x_min, x_max = axis.get_xlim()
    y_min, y_max = axis.get_ylim()
    x_span = (x_max - x_min) or 1.0
    y_span = (y_max - y_min) or 1.0
    placed: list[tuple] = []
    for _, row in labelled.nsmallest(label_top * 4, "effect_difference").iterrows():
        if len(placed) >= label_top:
            break
        fraction_x = (row["effect_difference"] - x_min) / x_span
        fraction_y = (row["neg_log10_p"] - y_min) / y_span
        if any(abs(fraction_x - px) < 0.05 and abs(fraction_y - py) < 0.028 for px, py in placed):
            continue
        placed.append((fraction_x, fraction_y))
        axis.annotate(
            row["gene"], (row["effect_difference"], row["neg_log10_p"]),
            textcoords="offset points", xytext=(5, 4), fontsize=9,
            color="#2e294e" if row["synlethdb_known"] else "0.15",
        )

    axis.set_xlabel(
        "Effect difference (altered mean - control mean gene effect)\n"
        "left = stronger dependency in the altered group"
    )
    axis.set_ylabel("$-\\log_{10}$(one-sided p)")
    axis.set_title(title, fontsize=12)
    axis.legend(loc="upper right", frameon=True, fontsize=9)
    figure.tight_layout()
    os.makedirs(os.path.dirname(os.path.abspath(output_path)), exist_ok=True)
    figure.savefig(output_path, dpi=200, bbox_inches="tight")
    plt.close(figure)
    return {
        "n_total": len(frame),
        "n_selected": int(selected.sum()),
        "n_nominal": int(nominal.sum()),
        "n_pan_essential": int(is_pan.sum()),
    }


def plot_dependency_boxplot(
    target_mutation: str,
    genes,
    output_path: str,
    histology: str | None = None,
    stratify_by: str = "mutation",
    data_lake_path: str | None = None,
    mutation_csv_path: str | None = None,
) -> str:
    """Plot per-cell-line CRISPR dependency of target genes in altered vs control ovarian lines.

    One box-and-strip panel per gene, with every cell line shown as a point so that a "difference"
    driven by one or two outlier lines is visible rather than hidden behind a mean. The groups are
    built exactly as in ``stratify_ovarian_cancer_dependency_by_mutation``, so a figure always
    matches the statistics reported there.

    Parameters
    ----------
    target_mutation : str
        HUGO symbol of the gene whose alteration defines the two groups, e.g. "BRCA1", "CCNE1".
    genes : list[str] | str
        Target genes to plot (list, or comma/whitespace separated string), e.g. the CANDIDATE_GENES
        line of a discovery run. Genes absent from the CRISPR library are reported and skipped.
    output_path : str
        Where to write the PNG.
    histology : str, optional
        Ovarian histology filter, e.g. "HGSOC", "clear cell". Default: all ovarian tumour lines.
    stratify_by : str, optional
        "mutation" (default), "amplification" or "deletion".
    data_lake_path : str, optional
        Directory holding the DepMap files (default: "./data/biomni_data/data_lake").
    mutation_csv_path : str, optional
        Genotype table overriding the default call source (see the analysis tool for the columns).

    Returns
    -------
    str
        A short log with the figure path, group sizes, genotype-call provenance and any skipped gene.

    """
    driver = str(target_mutation).strip().upper()
    mode = _STRATIFY_ALIASES.get(str(stratify_by).strip().lower())
    if mode is None:
        return f"FAILURE: stratify_by='{stratify_by}' is not recognised. Use one of {sorted(set(_STRATIFY_ALIASES.values()))}."
    altered_label, control_label, event_name, _ = _STRATIFY_MODES[mode]

    gene_list = _parse_gene_list(genes)
    if not gene_list:
        return "FAILURE: no gene was requested. Pass `genes` as a list or a comma separated string."

    try:
        groups = _ovarian_groups(driver, histology, mode, data_lake_path, mutation_csv_path)
    except (FileNotFoundError, ValueError, RuntimeError) as e:
        return f"FAILURE: {e}"

    mutant_ids, wildtype_ids = groups["mutant_ids"], groups["wildtype_ids"]
    if len(mutant_ids) < 2 or len(wildtype_ids) < 2:
        return (
            f"FAILURE: insufficient group sizes ({altered_label.lower()} n={len(mutant_ids)}, "
            f"{control_label.lower()} n={len(wildtype_ids)}) to plot a comparison."
        )

    cohort_title = "ovarian" + (f" ({histology})" if histology else "")
    try:
        plotted, skipped = _draw_dependency_boxplot(
            groups["bundle"], groups["cohort"], mutant_ids, wildtype_ids, gene_list,
            (altered_label, control_label), output_path,
            f"CRISPR dependency by {driver} {event_name} status - {cohort_title} "
            f"({altered_label} n={len(mutant_ids)} vs {control_label} n={len(wildtype_ids)})",
        )
    except ValueError as e:
        return f"FAILURE: {e}"

    log = [
        f"Boxplot written to {output_path}",
        f"  Genes plotted: {', '.join(plotted)}",
        f"  Groups: {altered_label} n={len(mutant_ids)} vs {control_label} n={len(wildtype_ids)}",
        f"  {event_name.capitalize()} calls: {groups['genotype_source']}",
    ]
    if skipped:
        log.append(f"  Skipped (not in the CRISPR library): {', '.join(skipped)}")
    log.append(
        "  Each point is one cell line; per-panel delta and one-sided Welch p are recomputed from the "
        "plotted values, so the figure and the numbers cannot drift apart."
    )
    return "\n".join(log)


def plot_dependency_volcano(
    results_csv_path: str,
    output_path: str,
    title: str | None = None,
    p_threshold: float = 0.05,
    fdr_threshold: float = 0.25,
    min_effect_difference: float = -0.2,
    exclude_pan_essential: bool = True,
    label_top: int = 15,
) -> str:
    """Draw a volcano plot from a dependency table written by the ovarian stratification tool.

    x is the effect difference (altered mean - control mean gene effect, so negative = the altered
    group is more dependent) and y is -log10 of the one-sided p-value. Genes passing the FDR gate
    are highlighted and labelled, genes that are only nominally significant are shown in a second
    colour, pan-essential genes are marked separately so that a common-essential gene is never
    mistaken for a selective vulnerability, and known SynLethDB partners are ringed.

    Parameters
    ----------
    results_csv_path : str
        CSV written by ``stratify_ovarian_cancer_dependency_by_mutation(output_csv_path=...)``, or
        any table with at least gene, effect_difference and p_one_sided columns.
    output_path : str
        Where to write the PNG.
    title : str, optional
        Figure title. Defaults to the CSV file name.
    p_threshold, fdr_threshold, min_effect_difference, exclude_pan_essential : optional
        Thresholds used to colour the points; pass the same values as the analysis run.
    label_top : int, optional
        Number of gene labels to draw, taken from the strongest effect differences (default: 15).

    Returns
    -------
    str
        A short log with the figure path and the point counts per category.

    """
    import pandas as pd

    if not os.path.exists(results_csv_path):
        return f"FAILURE: {results_csv_path} does not exist. Run the analysis with output_csv_path first."
    frame = pd.read_csv(results_csv_path)
    try:
        counts = _draw_volcano(
            frame, output_path, title or f"Differential dependency - {os.path.basename(results_csv_path)}",
            p_threshold, fdr_threshold, min_effect_difference, exclude_pan_essential, label_top,
        )
    except ValueError as e:
        return f"FAILURE: {e}"

    return "\n".join(
        [
            f"Volcano plot written to {output_path}",
            f"  Genes plotted: {counts['n_total']}",
            f"  Significant (q < {fdr_threshold}): {counts['n_selected']}",
            f"  Nominal only (p < {p_threshold}): {counts['n_nominal']}",
            f"  Pan-essential, marked separately: {counts['n_pan_essential']}",
            f"  Source table: {results_csv_path}",
        ]
    )
