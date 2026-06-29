"""Callbacks for contextual controls in sidebar (tab-specific options)."""

from __future__ import annotations

import logging
from typing import Dict

import dash_bootstrap_components as dbc
from dash import Input, Output, State, dcc, html

from .apis import (
    get_bulk_multi_api,
    get_expression_api,
    get_gene_search_api,
    get_landscape_view_api,
    get_scatac_api,
    get_scrna_api,
    is_bulkmulti_data_available,
)
from .species_data import bulkmulti_data_missing_alert, get_available_datasets


def register(app, *, logger: logging.Logger) -> None:
    """Register contextual controls callbacks."""

    @app.callback(
        Output("contextual-controls-container", "children"),
        Input("main-tabs", "value"),
        State("current-species-store", "data")
    )
    def update_contextual_controls(tab, species_key):
        logger.info(f"🔄 update_contextual_controls called with: {tab}")
        if tab == "tab-summary":
            # Get column configuration from adapter instead of hardcoded values
            current_gene_search_api = get_gene_search_api(species_key)
            if current_gene_search_api:
                columns_config = current_gene_search_api.adapter.get_summary_columns()

                # Use all available columns as default, ordered by the config
                default_columns = list(columns_config.keys())

                # Create options in the order defined by the adapter
                ordered_options = [
                    {'label': columns_config[col], 'value': col}
                    for col in default_columns
                ]
            else:
                # Generic fallback if no API available
                columns_config = {
                    'gene_name': 'Gene Name',
                    'primary_id': 'Gene ID',
                    'locus': 'Locus',
                    'tf_families': 'TF Family',
                    'description': 'Description'
                }
                default_columns = list(columns_config.keys())
                ordered_options = [
                    {'label': columns_config[col], 'value': col}
                    for col in default_columns
                ]

            return html.Div([
                html.Label("Summary Table Columns", className="form-label fw-bold"),
                dbc.Checklist(options=ordered_options, value=default_columns, id="summary-column-checklist"),
            ])
        elif tab == "tab-bulk-rna":
            current_expression_api = get_expression_api(species_key)
            if current_expression_api is None:
                return dbc.Alert("No bulk RNA expression data available for this species.", color="warning", className="mt-2")

            # Get available sources for selection
            sources = current_expression_api.get_available_sources()
            source_options = [{'label': f"{s['source']} ({s['sample_count']} samples)", 'value': s['source']} for s in sources]

            return html.Div([
                html.Label("Study", className="form-label fw-bold mb-2"),
                dbc.Checklist(
                    id="expression-sources-checklist",
                    options=source_options,
                    value=[sources[0]['source']] if sources else [],  # Default to first source
                    className="mb-3"
                ),
                html.Label("Sample Selection", className="form-label fw-bold mb-2"),
                html.Div([
                    dbc.Switch(
                        id="custom-sample-selection-switch",
                        label="Custom sample selection",
                        value=False,
                        className="mb-2"
                    ),
                    html.Div(id="sample-selection-container", className="mb-3")
                ]),
                html.Label("Visualization Type", className="form-label fw-bold mb-2"),
                dbc.RadioItems(
                    id="expression-viz-type-radio",
                    options=[
                        {'label': 'Heatmap', 'value': 'heatmap'},
                        {'label': 'Line Plot', 'value': 'lineplot'},
                    ],
                    value='heatmap',
                    className="mb-3"
                ),
                html.Label("Data Transformation", className="form-label fw-bold mb-2"),
                dbc.RadioItems(
                    id="expression-transform-radio",
                    options=[
                        {'label': 'Raw TPM', 'value': 'raw'},
                        {'label': 'Log2(TPM + 1)', 'value': 'log2'},
                        {'label': 'Z-score (by gene)', 'value': 'zscore_rows'},
                    ],
                    value='log2',
                    className="mb-3"
                ),
                html.Label("Display Options", className="form-label fw-bold mb-2"),
                dbc.Switch(
                    id="expression-cluster-switch",
                    label="Cluster genes",
                    value=False,
                    className="mb-2"
                ),
                dbc.Switch(
                    id="expression-show-values-switch",
                    label="Show values on heatmap",
                    value=False,
                    className="mb-2"
                ),
                html.Label("Export Data", className="form-label fw-bold mb-2 mt-4"),
                dbc.Button("Download CSV", id="export-csv-btn", color="secondary", size="sm", className="w-100"),
                dcc.Download(id="download-expression-data")
            ])
        elif tab == "tab-scrna":
            current_scrna_api = get_scrna_api(species_key)
            if current_scrna_api is None:
                return dbc.Alert("No single-cell RNA-seq data available for this species.", color="warning", className="mt-2")

            # Get available studies
            studies = current_scrna_api.get_available_studies()
            study_options = [{'label': study['name'], 'value': study['id']} for study in studies]
            default_study = studies[0]['id'] if studies else None

            return html.Div([
                # Study selection
                html.Label("Study", className="form-label fw-bold mb-2"),
                dbc.Select(
                    id="scrna-study-select",
                    options=study_options,
                    value=default_study,
                    placeholder="Select a study...",
                    className="mb-3"
                ),

                # Sample selection (will be populated by callback)
                html.Label("Sample", className="form-label fw-bold mb-2"),
                dbc.Select(
                    id="scrna-sample-select",
                    placeholder="Select a sample...",
                    className="mb-3"
                ),

                # Visualization type
                html.Label("Visualization Type", className="form-label fw-bold mb-2"),
                dbc.RadioItems(
                    id="scrna-viz-type-radio",
                    options=[
                        {'label': 'Expression UMAP', 'value': 'umap'},
                        {'label': 'Cluster UMAP', 'value': 'cluster'},
                        {'label': 'Cluster Violin', 'value': 'violin'},
                    ],
                    value='umap',
                    className="mb-3"
                ),

                # Data version selection
                html.Label("Data Version", className="form-label fw-bold mb-2"),
                dbc.RadioItems(
                    id="scrna-data-version-radio",
                    options=[
                        {'label': '📊 Normalized log1p', 'value': 'norm'},
                        {'label': '🪄 MAGIC Smoothed', 'value': 'magic'}
                    ],
                    value='norm',
                    className="mb-3"
                ),

                # Clustering method (will be populated by callback)
                html.Label("Clustering Method", className="form-label fw-bold mb-2"),
                dbc.Select(
                    id="scrna-cluster-select",
                    placeholder="Select clustering method...",
                    className="mb-3",
                    disabled=True  # Default disabled until a cluster-based viz is chosen
                ),

                # Color scale
                html.Label("Color Scale", className="form-label fw-bold mb-2"),
                dbc.Select(
                    id="scrna-colorscale-select",
                    options=[
                        {'label': 'Viridis', 'value': 'viridis'},
                        {'label': 'YlGnBu', 'value': 'YlGnBu'},
                        {'label': 'Plasma', 'value': 'plasma'},
                        {'label': 'Blues', 'value': 'blues'},
                        {'label': 'Oranges', 'value': 'oranges'},
                        {'label': 'Reds', 'value': 'reds'}
                    ],
                    value='oranges',
                    className="mb-3"
                ),

                # Expression threshold
                html.Label("Expression Threshold", className="form-label fw-bold mb-2 mt-4"),
                dcc.Slider(
                    id="scrna-threshold-slider",
                    min=0,
                    max=5,
                    step=0.1,
                    value=0.1,
                    marks={i: str(i) for i in range(6)},
                    tooltip={"placement": "bottom", "always_visible": True},
                    className="mb-3"
                )
            ])
        elif tab == "tab-scatac":
            current_scatac_api = get_scatac_api(species_key)
            if current_scatac_api is None:
                return dbc.Alert("No single-cell ATAC-seq data available for this species.", color="warning", className="mt-2")

            # Get available studies
            studies = current_scatac_api.get_available_studies()
            study_options = [{'label': study['name'], 'value': study['id']} for study in studies]
            default_study = studies[0]['id'] if studies else None

            return html.Div([
                # Study selection
                html.Label("Study", className="form-label fw-bold mb-2"),
                dbc.Select(
                    id="scatac-study-select",
                    options=study_options,
                    value=default_study,
                    placeholder="Select a study...",
                    className="mb-3"
                ),

                # Sample selection (will be populated by callback)
                html.Label("Sample", className="form-label fw-bold mb-2"),
                dbc.Select(
                    id="scatac-sample-select",
                    placeholder="Select a sample...",
                    className="mb-3"
                ),

                # Visualization type
                html.Label("Visualization Type", className="form-label fw-bold mb-2"),
                dbc.RadioItems(
                    id="scatac-viz-type-radio",
                    options=[
                        {'label': 'Gene Activity UMAP', 'value': 'umap'},
                        {'label': 'Cluster UMAP', 'value': 'cluster'},
                        {'label': 'Gene Activity Violin', 'value': 'stats'},
                    ],
                    value='umap',
                    className="mb-3"
                ),

                # Data version (placeholder for future multi-version GAS)
                html.Label("Data Version", className="form-label fw-bold mb-2"),
                dbc.RadioItems(
                    id="scatac-data-version-radio",
                    options=[
                        {'label': '📊 Gene Activity (log1p)', 'value': 'activity'}
                    ],
                    value='activity',
                    className="mb-3"
                ),

                # Clustering method (will be populated by callback)
                html.Label("Clustering Method", className="form-label fw-bold mb-2"),
                dbc.Select(
                    id="scatac-cluster-select",
                    placeholder="Select clustering method...",
                    className="mb-3",
                    disabled=True  # Default disabled until a cluster-based viz is chosen
                ),

                # Color scale
                html.Label("Color Scale", className="form-label fw-bold mb-2"),
                dbc.Select(
                    id="scatac-colorscale-select",
                    options=[
                        {'label': 'Viridis', 'value': 'viridis'},
                        {'label': 'YlGnBu', 'value': 'YlGnBu'},
                        {'label': 'Plasma', 'value': 'plasma'},
                        {'label': 'Blues', 'value': 'blues'},
                        {'label': 'Oranges', 'value': 'oranges'},
                        {'label': 'Reds', 'value': 'reds'}
                    ],
                    value='oranges',
                    className="mb-3"
                ),
                html.Label("Gene Activity Threshold", className="form-label fw-bold mb-2 mt-4"),
                dcc.Slider(
                    id="scatac-threshold-slider",
                    min=0,
                    max=2,
                    step=0.1,
                    value=0.1,
                    marks={i*0.5: str(i*0.5) for i in range(5)},
                    tooltip={"placement": "bottom", "always_visible": True},
                    className="mb-3"
                )
            ])
        elif tab == "tab-bulkmulti":
            bulk_multi_api = get_bulk_multi_api(species_key)
            logger.info(f"🔧 Creating BulkMulti controls, bulk_multi_api: {bulk_multi_api}")
            if bulk_multi_api is None:
                logger.error("❌ BulkMulti API is None!")
                return dbc.Alert("No multi-omics data available for this species.", color="warning", className="mt-2")
            if not is_bulkmulti_data_available(species_key):
                logger.error("❌ BulkMulti dataset missing on disk.")
                return bulkmulti_data_missing_alert()


            # Get available datasets
            if not species_key:
                logger.warning("⚠️ BulkMulti controls requested without active species context")
                return dbc.Alert("Select a species to configure BulkMulti options.", color="info", className="mt-2")
            available_datasets = get_available_datasets(species_key)
            dataset_options = [{'label': dataset, 'value': dataset} for dataset in available_datasets]
            default_dataset = available_datasets[0] if available_datasets else None
            if not dataset_options:
                return dbc.Alert("No BulkMulti datasets detected under the current data directory.", color="warning", className="mt-2")

            return html.Div([
                # 1. Study Selection
                html.Label("Study", className="form-label fw-bold mb-2"),
                dbc.Select(
                    id="bulkmulti-dataset-select",
                    options=dataset_options,
                    value=default_dataset,
                    className="mb-3"
                ),

                # 2. View Type
                html.Label("View Type", className="form-label fw-bold mb-2 mt-4"),
                dbc.RadioItems(
                    id="p2g-view-type-radio",
                    options=[
                        {'label': 'Linkage Statistics', 'value': 'statistics'},
                        {'label': 'Landscape View', 'value': 'landscape'},
                    ],
                    value='statistics',
                    className="mb-3"
                ),

                # 3. Dynamic view-specific options container
                html.Div(id="bulkmulti-view-options-container"),

                # 4. Statistics summary
            ])
        return dbc.Alert([
            dcc.Markdown("📋 **No additional options** available for this tab.", className="markdown-container mb-0")
        ], color="info", className="d-flex align-items-center")

    @app.callback(
        Output("sample-selection-container", "children"),
        [Input("custom-sample-selection-switch", "value"),
         Input("expression-sources-checklist", "value")],
        State("current-species-store", "data")
    )
    def update_sample_selection_container(custom_selection, selected_sources, species_key):
        current_expression_api = get_expression_api(species_key)
        if not custom_selection or not selected_sources or current_expression_api is None:
            return html.P("Using all samples from selected sources", className="text-muted small")

        # Get samples for selected sources
        samples = current_expression_api.get_samples_by_sources(selected_sources)
        if not samples:
            return dbc.Alert("No samples found for selected sources", color="warning", className="small")

        # Group samples by source
        samples_by_source: Dict[str, list] = {}
        for sample in samples:
            source = sample['source']
            samples_by_source.setdefault(source, []).append(sample)

        # Create expandable sections for each source
        accordion_items = []
        for source, source_samples in samples_by_source.items():
            sample_options = [
                {'label': s['sample'],
                 'value': f"{source}::{s['sample']}"}
                for s in source_samples
            ]

            accordion_items.append(
                dbc.AccordionItem([
                    dbc.Checklist(
                        id={"type": "sample-checklist", "source": source},
                        options=sample_options,
                        value=[f"{source}::{s['sample']}" for s in source_samples],  # Default: select all uniquely
                        className="small"
                    )
                ], title=f"{source} ({len(source_samples)} samples)")
            )

        return dbc.Accordion(accordion_items, start_collapsed=True, className="small")
