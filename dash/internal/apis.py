from __future__ import annotations

import logging
from pathlib import Path
from threading import Lock
from typing import Any, Dict, Optional

from backend.bulk_multi_api import BulkMultiAPI
from backend.expression_api import ExpressionAPI
from backend.gene_search import GeneSearchAPI
from backend.landscape_view_api import LandscapeViewAPI
from backend.scatac_api import scATACAPI
from backend.scrna_api import scRNAAPI

from .species_data import validate_species_data


logger = logging.getLogger("miniodp.dash")


# Global variables for species-specific APIs
species_api_cache: Dict[tuple[str, Optional[str]], Dict[str, Any]] = {}
species_cache_lock = Lock()


def _normalize_dataset_key(dataset: Optional[str]) -> Optional[str]:
    if dataset is None:
        return None
    if not isinstance(dataset, str):
        raise TypeError("Dataset must be a string or None")
    normalized = dataset.strip()
    return normalized or None


def initialize_apis(species_key: str, dataset: Optional[str] = None) -> Optional[Dict[str, Any]]:
    """Initialize all API instances for a specific species."""
    if not species_key:
        raise ValueError("Species key is required and cannot be empty")
    if not isinstance(species_key, str):
        raise TypeError("Species key must be a string")

    species_key = species_key.strip()
    if not species_key:
        raise ValueError("Species key cannot be empty or whitespace")

    # Validate species data first
    validation = validate_species_data(species_key)
    if not validation["valid"]:
        logger.error(f"❌ Species {species_key} validation failed:")
        for missing in validation["missing_required"]:
            logger.info(f"  - Missing: {missing}")
        return None

    logger.info(f"📋 Species {species_key} validation passed")
    logger.info(f"  Available features: {', '.join(validation['available_features'])}")
    if validation["missing_optional"]:
        logger.info(f"  Optional features not available: {', '.join(validation['missing_optional'])}")

    normalized_dataset = _normalize_dataset_key(dataset)

    try:
        apis: Dict[str, Any] = {}

        # Create species adapter that will be shared across APIs
        from backend.species_adapters import create_species_adapter

        adapter = create_species_adapter(species_key)

        # Always initialize core APIs with adapter
        apis["gene_search_api"] = GeneSearchAPI(species_key=species_key)

        # Initialize optional APIs based on available data
        if "Bulk RNA-seq data" in validation["available_features"]:
            apis["expression_api"] = ExpressionAPI(species_key=species_key, adapter=adapter)

        if "Single-cell RNA data" in validation["available_features"]:
            apis["scrna_api"] = scRNAAPI(species_key=species_key)

        if "Single-cell ATAC data" in validation["available_features"]:
            apis["scatac_api"] = scATACAPI(species_key=species_key)

        if "Bulk multi-omics analysis data" in validation["available_features"]:
            apis["bulk_multi_api"] = BulkMultiAPI(species_key=species_key, dataset=normalized_dataset)
            apis["landscape_view_api"] = LandscapeViewAPI(species_key=species_key, dataset=normalized_dataset)

        logger.info(f"✅ APIs initialized successfully for species: {species_key}")
        logger.info(f"  Loaded APIs: {', '.join(apis.keys())}")
        return apis
    except Exception as exc:
        logger.error(f"❌ Failed to initialize APIs for {species_key}: {exc}")
        return None


def _normalize_species_key(species_key: Optional[str]) -> Optional[str]:
    if not species_key or not isinstance(species_key, str):
        return None
    normalized = species_key.strip()
    return normalized or None


def get_species_apis(
    species_key: Optional[str],
    dataset: Optional[str] = None,
) -> Optional[Dict[str, Any]]:
    """Return cached API instances for a species, initializing if needed."""
    normalized = _normalize_species_key(species_key)
    if not normalized:
        return None
    normalized_dataset = _normalize_dataset_key(dataset)
    cache_key = (normalized, normalized_dataset)
    apis = species_api_cache.get(cache_key)
    if apis:
        return apis
    with species_cache_lock:
        apis = species_api_cache.get(cache_key)
        if apis:
            return apis
        apis = initialize_apis(normalized, normalized_dataset)
        if apis:
            species_api_cache[cache_key] = apis
            return apis
    return None


def get_api_for_species(species_key: Optional[str], api_name: str, dataset: Optional[str] = None):
    normalized = _normalize_species_key(species_key)
    if not normalized:
        return None

    normalized_dataset = _normalize_dataset_key(dataset)
    cache_key = (normalized, normalized_dataset)

    apis = get_species_apis(normalized, normalized_dataset)
    if not apis:
        return None

    api = apis.get(api_name)
    if api is not None:
        return api

    logger.warning(
        "API '%s' missing from cache for species '%s' (dataset=%s); rebuilding cache entry",
        api_name,
        normalized,
        normalized_dataset,
    )

    with species_cache_lock:
        species_api_cache.pop(cache_key, None)
        refreshed = initialize_apis(normalized, normalized_dataset)
        if refreshed:
            species_api_cache[cache_key] = refreshed
            return refreshed.get(api_name)

    return None


def get_gene_search_api(species_key: Optional[str], dataset: Optional[str] = None):
    return get_api_for_species(species_key, "gene_search_api", dataset)


def get_expression_api(species_key: Optional[str], dataset: Optional[str] = None):
    return get_api_for_species(species_key, "expression_api", dataset)


def get_scrna_api(species_key: Optional[str], dataset: Optional[str] = None):
    return get_api_for_species(species_key, "scrna_api", dataset)


def get_scatac_api(species_key: Optional[str], dataset: Optional[str] = None):
    return get_api_for_species(species_key, "scatac_api", dataset)


def get_bulk_multi_api(species_key: Optional[str], dataset: Optional[str] = None):
    return get_api_for_species(species_key, "bulk_multi_api", dataset)


def get_landscape_view_api(species_key: Optional[str], dataset: Optional[str] = None):
    return get_api_for_species(species_key, "landscape_view_api", dataset)


def is_bulkmulti_data_available(species_key: Optional[str], dataset: Optional[str] = None) -> bool:
    """Check if BulkMulti API has a valid dataset on disk."""
    bulk_api = get_bulk_multi_api(species_key, dataset)
    if not bulk_api:
        return False
    db_path = getattr(bulk_api, "db_path", None)
    if not db_path:
        return False
    return Path(db_path).exists()
