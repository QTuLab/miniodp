from pathlib import Path
from typing import Any, Dict

import dash_bootstrap_components as dbc
from dash import html

from .settings import DASH_DATA_PATH


def get_available_species() -> list[str]:
    """Get list of available species from data directory."""
    data_path = Path(DASH_DATA_PATH)
    if not data_path.exists():
        raise FileNotFoundError(f"Data directory not found: {data_path}")

    species: list[str] = []
    for item in data_path.iterdir():
        if item.is_dir() and item.name != "common":
            species.append(item.name)

    if not species:
        raise ValueError(f"No species data found in directory: {data_path}")

    return sorted(species)


def validate_species_data(species_key: str) -> Dict[str, Any]:
    """Validate that required data exists for a species."""
    base_dir = Path(DASH_DATA_PATH) / species_key

    required_files = {
        "geneinfo.db": "Gene information database",
        "geneinfo.toml": "Gene search configuration",
    }

    optional_dirs = {
        "bulkRNA": "Bulk RNA-seq data",
        "scRNA": "Single-cell RNA data",
        "scATAC": "Single-cell ATAC data",
        "BulkMulti": "Bulk multi-omics analysis data",
    }

    validation_results: Dict[str, Any] = {
        "species_key": species_key,
        "valid": True,
        "missing_required": [],
        "missing_optional": [],
        "available_features": [],
    }

    # Check required files
    for file_name, description in required_files.items():
        file_path = base_dir / file_name
        if not file_path.exists():
            validation_results["missing_required"].append(f"{description} ({file_name})")
            validation_results["valid"] = False

    # Check optional directories
    for dir_name, description in optional_dirs.items():
        dir_path = base_dir / dir_name
        if dir_path.exists() and dir_path.is_dir():
            validation_results["available_features"].append(description)
        else:
            validation_results["missing_optional"].append(f"{description} ({dir_name})")

    return validation_results


def get_available_datasets(species_key: str) -> list[str]:
    """
    Discover available BulkMulti datasets for a given species.
    Returns list of dataset names sorted in reverse order (newest first).
    """
    if not species_key:
        raise ValueError("Species key is required and cannot be empty")
    if not isinstance(species_key, str):
        raise TypeError("Species key must be a string")

    species_key = species_key.strip()
    if not species_key:
        raise ValueError("Species key cannot be empty or whitespace")

    bulkmulti_path = Path(DASH_DATA_PATH) / species_key / "BulkMulti"
    if not bulkmulti_path.exists():
        return []

    datasets: list[str] = []
    for item in bulkmulti_path.iterdir():
        if item.is_dir():
            datasets.append(item.name)

    # Sort in reverse order (newest first, e.g., 2025Bulk before 2024Bulk)
    return sorted(datasets, reverse=True)


def bulkmulti_data_missing_alert():
    data_root = Path(DASH_DATA_PATH).resolve()
    return dbc.Alert(
        [
            html.H5("BulkMulti data is not available yet", className="alert-heading mb-1"),
            html.P(
                [
                    "Confirm that the data directory contains ",
                    html.Code("Species/BulkMulti/<dataset>/BulkMulti.db"),
                    " and the required BigWig/peaks files. Current search path:",
                    html.Br(),
                    html.Code(str(data_root)),
                ],
                className="mb-2",
            ),
            html.P("Mount or copy the dataset, then refresh the page.", className="mb-0"),
        ],
        color="warning",
        className="mt-3",
    )


def get_summary_columns(species_key: str, logger):
    """Get summary columns dynamically from species adapter."""
    try:
        from backend.species_adapters.factory import create_species_adapter

        adapter = create_species_adapter(species_key)
        columns = adapter.get_summary_columns()
        return list(columns.keys())
    except Exception as exc:
        logger.error(f"Error getting summary columns for {species_key}: {exc}")
        # Fallback to generic columns
        return ["gene_name", "primary_id", "locus", "tf_families", "description"]


def get_column_display_config(species_key: str, logger):
    """Get column display configuration from species adapter."""
    try:
        from backend.species_adapters.factory import create_species_adapter

        adapter = create_species_adapter(species_key)
        return adapter.get_column_display_config()
    except Exception as exc:
        logger.error(f"Error getting column display config for {species_key}: {exc}")
        return {}

