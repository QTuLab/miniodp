from __future__ import annotations

import logging
from typing import Sequence

import dash_bootstrap_components as dbc
from dash import Input, Output, html

from .apis import get_gene_search_api, get_species_apis
from .layout import (
    create_main_content,
    create_sidebar,
    create_species_error_layout,
    create_species_selection_layout,
)
from .settings import strip_dash_prefix


def _get_species_from_pathname(pathname: str | None) -> str | None:
    """Extract species key from URL pathname.

    With url_base_pathname configured, Dash may already strip the prefix so pathname can be:
    - / (root)
    - /zebrafish (species page)
    """
    if not pathname:
        return None

    normalized_path = strip_dash_prefix(pathname)
    if normalized_path in ("", "/"):
        return None

    path_parts = [p for p in normalized_path.strip("/").split("/") if p]
    return path_parts[0] if path_parts else None


def register(
    app,
    *,
    logger: logging.Logger,
    available_species: Sequence[str],
    tf_families_list: str,
    default_search_field_options,
) -> None:
    """Register routing-related callbacks."""

    @app.callback(
        [Output("page-content", "children"), Output("current-species-store", "data")],
        [Input("url", "pathname")],
    )
    def display_page(pathname):
        """Main routing callback to handle species-specific URLs."""
        species_key = _get_species_from_pathname(pathname)

        if species_key is None:
            return create_species_selection_layout(list(available_species)), None
        if species_key not in available_species:
            return create_species_error_layout(species_key, list(available_species)), species_key

        apis = get_species_apis(species_key)
        if not apis:
            return create_species_error_layout(species_key, list(available_species)), species_key

        logger.info(f"🔄 Routing to species: {species_key}")

        layout = html.Div(
            [create_sidebar(tf_families_list, default_search_field_options), create_main_content()]
        )
        return layout, species_key

    @app.callback(Output("current-species-display", "children"), Input("current-species-store", "data"))
    def update_species_display(species_key):
        if species_key:
            return dbc.Badge(
                f"Species: {species_key.title()}",
                color="primary",
                className="fs-5 fw-bold px-3 py-2",
            )
        return ""

    @app.callback(Output("search-field-select", "options"), Input("current-species-store", "data"))
    def update_search_field_options(species_key):
        current_gene_search_api = get_gene_search_api(species_key)
        if species_key and current_gene_search_api:
            return current_gene_search_api.get_available_search_fields()
        return default_search_field_options

