"""Callbacks for BulkMulti (multi-omics) analysis."""

from __future__ import annotations

import logging
import time

import dash
import dash_bootstrap_components as dbc
from dash import Input, Output, State, callback_context, dcc, html
from dash.exceptions import PreventUpdate

from .apis import (
    get_bulk_multi_api,
    get_landscape_view_api,
    is_bulkmulti_data_available,
)
from .plots import (
    _build_bulkmulti_view,
    convert_gene_ids_to_names,
    create_bulkmulti_landscape_view,
    create_bulkmulti_statistics_view,
    create_single_gene_landscape_view,
)
from .species_data import bulkmulti_data_missing_alert


def register(app, *, logger: logging.Logger) -> None:
    """Register BulkMulti callbacks."""

    # --- BulkMulti Dynamic Options Callback ---

    @app.callback(
        Output("bulkmulti-view-options-container", "children"),
        Input("p2g-view-type-radio", "value"),
        prevent_initial_call=False
    )
    def update_bulkmulti_view_options(view_type):
        """Update view-specific options based on selected view type"""
        if view_type == "landscape":
            return html.Div([
                html.Label("Landscape View Options", className="form-label fw-bold mb-2 text-primary"),

                # Region extension settings
                html.Label("Region Extension", className="form-label fw-bold mb-2"),
                dbc.Row([
                    dbc.Col([
                        html.Label("Upstream (kb)", className="form-label small"),
                        dcc.Input(
                            id="landscape-upstream-input",
                            type="number",
                            value=50,
                            min=10,
                            max=500,
                            step=10,
                            className="form-control form-control-sm"
                        )
                    ], width=6),
                    dbc.Col([
                        html.Label("Downstream (kb)", className="form-label small"),
                        dcc.Input(
                            id="landscape-downstream-input",
                            type="number",
                            value=50,
                            min=10,
                            max=500,
                            step=10,
                            className="form-control form-control-sm"
                        )
                    ], width=6)
                ], className="mb-3"),

                # ELS guide lines option
                html.Label("Display Options", className="form-label fw-bold mb-2"),
                dbc.Switch(
                    id="landscape-els-lines-switch",
                    label="Show ELS guide lines on ATAC tracks",
                    value=False,
                    className="mb-3"
                ),


                # ATAC visualization settings
                html.Label("ATAC Visualization", className="form-label fw-bold mb-2"),
                html.Div([
                    html.Label("Signal Height Ratio", className="form-label small mb-2"),
                    html.P("Adjust the height ratio of ATAC signals relative to other tracks",
                           className="text-muted small mb-2"),
                    dcc.Slider(
                        id="landscape-atac-height-slider",
                        min=0.05,
                        max=1.0,  # Cap ratio at 1.0
                        step=0.05,
                        value=0.2,
                        marks={i*0.2: f"{i*0.2:.1f}" for i in range(6)},  # Align slider marks with the range
                        tooltip={"placement": "bottom", "always_visible": True},
                        className="mb-2"
                    )
                ], className="mb-3", style={
                    "backgroundColor": "#f8f9fa",  # Light gray background
                    "padding": "15px",
                    "borderRadius": "8px",
                    "border": "1px solid #e9ecef"
                }),

                # Linkage Strength for landscape view
                html.Label("Linkage Strength", className="form-label fw-bold mb-2"),
                dcc.Slider(
                    id="landscape-value-threshold-slider",
                    min=0,
                    max=1,
                    step=0.05,
                    value=0.2,
                    marks={i*0.2: f"{i*0.2:.1f}" for i in range(6)},
                    tooltip={"placement": "bottom", "always_visible": True},
                    className="mb-3"
                ),

                # Apply changes button
                html.Div([
                    dcc.Markdown("🔄 **Apply Changes** - Click to update the landscape view with your new settings:",
                           className="markdown-container small mb-2"),
                    dbc.Button(
                        "Apply Changes",
                        id="landscape-apply-changes-btn",
                        color="primary",
                        size="sm",
                        className="mt-1",
                        style={"width": "100%"}
                    )
                ], className="mt-3")
            ])

        elif view_type == "statistics":
            return html.Div([
                html.Label("Statistics Options", className="form-label fw-bold mb-2 text-success"),
                html.P("Statistical analysis options:", className="small text-muted mb-2"),

                # Linkage Strength for statistics view
                html.Label("Linkage Strength", className="form-label fw-bold mb-2"),
                dcc.Slider(
                    id="statistics-value-threshold-slider",
                    min=0,
                    max=1,
                    step=0.05,
                    value=0.2,
                    marks={i*0.2: f"{i*0.2:.1f}" for i in range(6)},
                    tooltip={"placement": "bottom", "always_visible": True},
                    className="mb-3"
                ),

                # Chart type selection
                html.Label("Chart Types", className="form-label fw-bold mb-2"),
                dbc.Checklist(
                    id="statistics-chart-types-checklist",
                    options=[
                        {'label': 'Linkage Count Bar Chart', 'value': 'bar_chart'},
                        {'label': 'Strength vs Distance Scatter', 'value': 'scatter_plot'},
                        {'label': 'Summary Statistics Table', 'value': 'stats_table'}
                    ],
                    value=['bar_chart', 'scatter_plot', 'stats_table'],
                    className="mb-3"
                )
            ])

        return html.Div([
            dbc.Alert([
                dcc.Markdown("🎯 **Select a view type** above to see available configuration options.", className="markdown-container mb-0")
            ], color="light", className="d-flex align-items-center")
        ])


    # --- BulkMulti Visualization Callback ---


    # --- Landscape View Options Update Callback ---
    @app.callback(
        Output("peak2gene-content", "children", allow_duplicate=True),
        Input("landscape-apply-changes-btn", "n_clicks"),
        [State("landscape-upstream-input", "value"),
         State("landscape-downstream-input", "value"),
         State("landscape-els-lines-switch", "value"),
         State("landscape-atac-height-slider", "value"),
         State("landscape-value-threshold-slider", "value"),
         State("target-gene-list-store", "data"),
         State("main-tabs", "value"),
         State("p2g-view-type-radio", "value"),
         State("bulkmulti-dataset-select", "value"),
         State("current-species-store", "data")],
        prevent_initial_call=True
    )
    def update_landscape_view_options(n_clicks, upstream_kb, downstream_kb, show_els_lines, atac_height_ratio,
                                    value_threshold, target_gene_ids, active_tab, view_type, selected_dataset, species_key):
        """Update landscape view when Apply Changes button is clicked"""
        if n_clicks is None:
            logger.info("🎛️ No button click detected")
            raise PreventUpdate

        logger.info(f"🎛️ Apply Changes clicked ({n_clicks} times):")
        logger.info(f"  - upstream_kb: {upstream_kb}")
        logger.info(f"  - downstream_kb: {downstream_kb}")
        logger.info(f"  - show_els_lines: {show_els_lines}")
        logger.info(f"  - atac_height_ratio: {atac_height_ratio}")

        bulk_multi_api = get_bulk_multi_api(species_key, selected_dataset)

        # Only update if we're on the BulkMulti tab and in landscape view
        if (active_tab != "tab-bulkmulti" or not target_gene_ids or not bulk_multi_api or
            view_type != "landscape"):
            logger.info("🎛️ Not updating - not in landscape view")
            raise PreventUpdate

        if not is_bulkmulti_data_available(species_key, selected_dataset):
            logger.error("❌ BulkMulti dataset missing when applying landscape changes.")
            return bulkmulti_data_missing_alert()

        try:
            gene_names = convert_gene_ids_to_names(species_key, target_gene_ids)
            if not gene_names:
                raise PreventUpdate

            # Create landscape config from user options
            landscape_config = {
                'upstream_kb': int(upstream_kb) if upstream_kb is not None else 50,
                'downstream_kb': int(downstream_kb) if downstream_kb is not None else 50,
                'show_els_guide_lines': bool(show_els_lines) if show_els_lines is not None else False,
                'atac_height_ratio': float(atac_height_ratio) if atac_height_ratio is not None else 0.2
            }

            # Safety bounds
            landscape_config['upstream_kb'] = max(10, min(500, landscape_config['upstream_kb']))
            landscape_config['downstream_kb'] = max(10, min(500, landscape_config['downstream_kb']))
            landscape_config['atac_height_ratio'] = max(0.05, min(1.0, landscape_config['atac_height_ratio']))

            logger.info(f"🎛️ Creating landscape view with config: {landscape_config}")
            return create_bulkmulti_landscape_view(bulk_multi_api, landscape_view_api, gene_names, value_threshold, species_key, landscape_config)

        except Exception as e:
            logger.error(f"❌ Error updating landscape options: {e}")
            raise PreventUpdate

    # --- Landscape View Pagination Callback ---
    @app.callback(
        [Output("peak2gene-content", "children", allow_duplicate=True),
         Output("landscape-page-store", "data", allow_duplicate=True)],
        [Input("landscape-prev-btn", "n_clicks"),
         Input("landscape-next-btn", "n_clicks")],
        [State("landscape-page-store", "data"),
         State("landscape-genes-store", "data"),
         State("landscape-config-store", "data"),
         State("landscape-threshold-store", "data"),
         State("bulkmulti-dataset-select", "value"),
         State("main-tabs", "value"),
         State("current-species-store", "data")],
        prevent_initial_call=True
    )
    def handle_landscape_pagination(prev_clicks, next_clicks, current_page, gene_names,
                                  landscape_config, value_threshold, selected_dataset, active_tab, species_key):
        """Handle landscape view pagination navigation"""
        logger.info(f"🔄 Pagination callback triggered:")
        logger.info(f"  - prev_clicks: {prev_clicks}")
        logger.info(f"  - next_clicks: {next_clicks}")
        logger.info(f"  - current_page: {current_page}")
        logger.info(f"  - gene_names: {gene_names}")
        logger.info(f"  - active_tab: {active_tab}")

        # Only handle if we're on BulkMulti tab
        bulk_multi_api = get_bulk_multi_api(species_key, selected_dataset)
        landscape_view_api = get_landscape_view_api(species_key, selected_dataset)
        if active_tab != "tab-bulkmulti":
            logger.error("❌ Not on BulkMulti tab, preventing update")
            raise PreventUpdate

        if not is_bulkmulti_data_available(species_key, selected_dataset):
            logger.error("❌ BulkMulti dataset missing during pagination.")
            return bulkmulti_data_missing_alert(), dash.no_update

        # Check if we have required data
        if not gene_names or current_page is None or not landscape_config:
            logger.error("❌ Missing required data for pagination")
            logger.info(f"  - gene_names valid: {bool(gene_names)}")
            logger.info(f"  - current_page valid: {current_page is not None}")
            logger.info(f"  - landscape_config valid: {bool(landscape_config)}")
            raise PreventUpdate

        ctx = callback_context
        if not ctx.triggered:
            logger.error("❌ No button triggered")
            raise PreventUpdate

        # Determine which button was clicked
        button_id = ctx.triggered[0]['prop_id'].split('.')[0]
        logger.info(f"🔄 Button clicked: {button_id}")

        new_page = current_page
        if button_id == "landscape-prev-btn" and current_page > 0:
            new_page = current_page - 1
            logger.info(f"🔄 Moving to previous page: {new_page}")
        elif button_id == "landscape-next-btn" and current_page < len(gene_names) - 1:
            new_page = current_page + 1
            logger.info(f"🔄 Moving to next page: {new_page}")
        else:
            logger.error(f"❌ Invalid navigation: page={current_page}, total={len(gene_names)}")
            raise PreventUpdate

        try:
            # Create single gene landscape view for the new page
            logger.info(f"🔄 Creating landscape view for gene {gene_names[new_page]} (page {new_page})")

            # Create complete paginated layout with updated stores
            content = html.Div([
                # Update all stores with current data
                dcc.Store(id="landscape-genes-store", data=gene_names),
                dcc.Store(id="landscape-config-store", data=landscape_config),
                dcc.Store(id="landscape-threshold-store", data=value_threshold),
                dcc.Store(id="landscape-page-store", data=new_page),

                # Create single gene view for the new page
                create_single_gene_landscape_view(
                    bulk_multi_api,
                    landscape_view_api,
                    gene_names[new_page],
                    new_page,
                    len(gene_names),
                    landscape_config,
                    value_threshold
                )
            ])

            return content, new_page

        except Exception as e:
            logger.error(f"❌ Error in pagination callback: {e}")
            import traceback
            traceback.print_exc()
            raise PreventUpdate

    @app.callback(
        [Output("peak2gene-content", "children"),
         Output("bulkmulti-render-request", "data")],
        [Input("target-gene-list-store", "data"),
         Input("main-tabs", "value"),
         Input("p2g-view-type-radio", "value"),
         Input("bulkmulti-dataset-select", "value")],
        # Statistics slider is handled by a dedicated callback
        [State("statistics-value-threshold-slider", "value"),
         State("current-species-store", "data")],
        prevent_initial_call=False
    )
    def update_peak2gene_content(target_gene_ids, active_tab, view_type, selected_dataset, value_threshold, species_key):
        """Prime BulkMulti render request and provide immediate loading feedback."""
        logger.info("🔍 update_peak2gene_content called (deferred mode)")
        logger.info(f"  - active_tab: {active_tab}")
        logger.info(f"  - target_gene_ids: {target_gene_ids}")
        logger.info(f"  - requested view_type: {view_type}")
        logger.info(f"  - selected_dataset: {selected_dataset}")
        logger.info(f"  - value_threshold slider: {value_threshold}")

        bulk_multi_api = get_bulk_multi_api(species_key, selected_dataset)
        landscape_view_api = get_landscape_view_api(species_key, selected_dataset)
        if active_tab != "tab-bulkmulti" or not target_gene_ids or not bulk_multi_api:
            logger.error("❌ BulkMulti tab inactive or missing prerequisites. Returning helpful message.")
            return (
                dbc.Alert([
                    dcc.Markdown("🧬 **Select target genes** to explore BulkMulti linkages and data visualization.", className="markdown-container mb-0")
                ], color="primary", className="d-flex align-items-center"),
                None
            )
        if not is_bulkmulti_data_available(species_key, selected_dataset):
            logger.error("❌ BulkMulti dataset missing, prompt user to upload data.")
            return bulkmulti_data_missing_alert(), None

        # Normalize defaults so downstream render logic receives consistent values
        if value_threshold is None:
            value_threshold = 0.2
            logger.info(f"🔧 Defaulting value_threshold to {value_threshold}")

        if not view_type:
            view_type = "statistics"
            logger.info(f"🔧 Defaulting view_type to {view_type}")

        # Provide immediate visual feedback while deferred rendering happens
        loading_view = html.Div([
            dbc.Spinner(color="primary", size="lg"),
            html.P("Loading BulkMulti view…", className="text-muted mt-3 mb-0")
        ], className="d-flex flex-column align-items-center justify-content-center py-5")

        request_payload = {
            "target_gene_ids": target_gene_ids,
            "view_type": view_type,
            "selected_dataset": selected_dataset,
            "value_threshold": value_threshold,
            "species_key": species_key,
            "ts": time.time()
        }

        logger.info(f"📨 Dispatching deferred render payload: {request_payload}")
        return loading_view, request_payload


    @app.callback(
        Output("peak2gene-content", "children", allow_duplicate=True),
        Input("bulkmulti-render-request", "data"),
        prevent_initial_call=True
    )
    def render_bulkmulti_content(request_payload):
        """Consume deferred BulkMulti render payload and produce final content."""
        if not request_payload:
            logger.info("ℹ️ No BulkMulti render payload received; preventing update.")
            raise PreventUpdate

        logger.info(f"📦 Processing deferred payload: {request_payload}")
        return _build_bulkmulti_view(
            request_payload.get("species_key"),
            request_payload.get("target_gene_ids", []),
            request_payload.get("view_type"),
            request_payload.get("selected_dataset"),
            request_payload.get("value_threshold", 0.2)
        )

    # --- Statistics View Real-time Update Callback ---
    @app.callback(
        Output("peak2gene-content", "children", allow_duplicate=True),
        Input("statistics-value-threshold-slider", "value"),
        [State("target-gene-list-store", "data"),
         State("main-tabs", "value"),
         State("p2g-view-type-radio", "value"),
         State("bulkmulti-dataset-select", "value"),
         State("current-species-store", "data")],
        prevent_initial_call=True
    )
    def update_statistics_view_realtime(value_threshold, target_gene_ids, active_tab, view_type, selected_dataset, species_key):
        """Real-time update for statistics view when threshold changes"""
        logger.info(f"🔍 Statistics real-time update: threshold={value_threshold}")

        # Only update if we're on BulkMulti tab and statistics view
        if active_tab != "tab-bulkmulti" or view_type != "statistics" or not target_gene_ids:
            raise PreventUpdate

        # Use default threshold if None
        if value_threshold is None:
            value_threshold = 0.2

        bulk_multi_api = get_bulk_multi_api(species_key, selected_dataset)
        if not bulk_multi_api:
            return dbc.Alert("BulkMulti API is not available for this species.", color="warning")

        try:
            # Convert gene names
            gene_names = convert_gene_ids_to_names(species_key, target_gene_ids)
            if not gene_names:
                raise PreventUpdate

            # Create statistics view
            stats = bulk_multi_api.get_linkage_statistics(gene_names)
            return create_bulkmulti_statistics_view(bulk_multi_api, gene_names, value_threshold, stats)

        except Exception as e:
            logger.error(f"❌ Error in statistics real-time update: {e}")
            raise PreventUpdate
