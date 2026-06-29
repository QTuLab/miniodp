"""Callbacks for single-cell ATAC-seq (scATAC) analysis."""

from __future__ import annotations

import logging
import traceback

import dash
import dash_bootstrap_components as dbc
from dash import Input, Output, State, callback_context, html, no_update

from .apis import get_scatac_api
from .plots import (
    create_accessibility_stats_plot,
    create_scatac_cluster_umap_plot,
    create_umap_accessibility_plot,
)


def register(app, *, logger: logging.Logger) -> None:
    """Register scATAC-seq callbacks."""

    # Study → Sample cascade callback
    @app.callback(
        [Output("scatac-sample-select", "options"),
         Output("scatac-sample-select", "value")],
        Input("scatac-study-select", "value"),
        State("current-species-store", "data")
    )
    def update_scatac_sample_options(study_id, species_key):
        """Update sample options based on selected study"""
        current_scatac_api = get_scatac_api(species_key)
        if not study_id or current_scatac_api is None:
            return [], None

        samples = current_scatac_api.get_samples_by_study(study_id)
        options = [{'label': f"{sample['display_name']} ({sample['metadata'].get('n_cells', 0):,} cells)",
                    'value': sample['id']} for sample in samples]
        default_value = samples[0]['id'] if samples else None

        return options, default_value

    # Sample → Cluster method cascade callback
    @app.callback(
        [Output("scatac-cluster-select", "options"),
         Output("scatac-cluster-select", "value")],
        Input("scatac-sample-select", "value"),
        State("current-species-store", "data")
    )
    def update_scatac_cluster_options(sample_id, species_key):
        """Update cluster method options based on selected sample"""
        current_scatac_api = get_scatac_api(species_key)
        if not sample_id or current_scatac_api is None:
            return [], None

        # Get cluster info from the selected sample
        cluster_info = current_scatac_api.get_cluster_info(sample_id)

        if cluster_info:
            options = [{'label': col.replace('_', ' ').title(), 'value': col} for col in cluster_info.keys()]

            def pick_preferred(opts):
                if not opts:
                    return None
                lowered = [(opt['value'], opt['label'].lower()) for opt in opts]
                for val, label in lowered:
                    if 'celltype' in label or 'cell type' in label or 'cell_type' in label:
                        return val
                for val, label in lowered:
                    if 'cluster' in label:
                        return val
                return opts[0]['value']

            default_value = pick_preferred(options)
        else:
            # Fallback options
            options = [
                {'label': 'ATAC Resolution 0.8', 'value': 'ATAC_snn_res.0.8'},
                {'label': 'Seurat Clusters', 'value': 'seurat_clusters'}
            ]
            default_value = 'ATAC_snn_res.0.8'

        return options, default_value

    # Toggle scATAC control states based on visualization type and sample readiness
    @app.callback(
        [Output("scatac-cluster-select", "disabled"),
         Output("scatac-colorscale-select", "disabled"),
         Output("scatac-threshold-slider", "disabled")],
        Input("scatac-viz-type-radio", "value"),
        Input("scatac-sample-select", "value")
    )
    def toggle_scatac_control_states(viz_type, sample_id):
        has_sample = bool(sample_id)
        cluster_disabled = (viz_type not in ("stats", "cluster")) or not has_sample
        color_disabled = (viz_type != "umap") or not has_sample
        threshold_disabled = (viz_type != "umap") or not has_sample
        return cluster_disabled, color_disabled, threshold_disabled

    @app.callback(
        Output("scatac-content", "children"),
        [Input("target-gene-list-store", "data"),
         Input("scatac-sample-select", "value"),
         Input("scatac-viz-type-radio", "value"),
         Input("scatac-data-version-radio", "value"),
         Input("scatac-cluster-select", "value"),
         Input("scatac-colorscale-select", "value"),
         Input("scatac-threshold-slider", "value"),
         Input("scatac-page-store", "data"),
         Input("main-tabs", "value")],
        State("current-species-store", "data")
    )
    def update_scatac_content(target_gene_ids, sample_id, viz_type, data_version, cluster_column, colorscale, threshold, current_page, active_tab, species_key):
        """Update scATAC-seq visualization content"""
        logger.info(f"🧬 update_scatac_content called:")
        logger.info(f"  - active_tab: {active_tab}")
        logger.info(f"  - target_gene_ids: {target_gene_ids}")
        logger.info(f"  - sample_id: {sample_id}")
        logger.info(f"  - viz_type: {viz_type}")
        logger.info(f"  - data_version: {data_version}")

        if active_tab != "tab-scatac":
            return no_update

        if not target_gene_ids:
            return dbc.Alert("Select genes from the sidebar to visualize chromatin accessibility.", color="info", className="mt-4")

        current_scatac_api = get_scatac_api(species_key)
        if current_scatac_api is None:
            return dbc.Alert("scATAC API is not available.", color="danger", className="mt-4")

        if not sample_id:
            return dbc.Alert("No sample selected.", color="warning", className="mt-4")

        try:
            genes_per_page = 4
            current_page = current_page or 0
            start_idx = current_page * genes_per_page
            end_idx = start_idx + genes_per_page
            current_genes = target_gene_ids[start_idx:end_idx]

            if not current_genes:
                return dbc.Alert("No genes to display on this page.", color="warning", className="mt-4")

            total_pages = (len(target_gene_ids) + genes_per_page - 1) // genes_per_page
            pagination_controls = html.Div([
                dbc.Row([
                    dbc.Col([
                        html.H6("scATAC-seq Visualization", className="mb-0")
                    ], width=6),
                    dbc.Col([
                        dbc.ButtonGroup([
                            dbc.Button("◀ Previous", id="scatac-prev-btn",
                                       disabled=current_page == 0, size="sm", color="secondary"),
                            dbc.Button(f"{current_page + 1} / {total_pages}",
                                       color="light", disabled=True, size="sm"),
                            dbc.Button("Next ▶", id="scatac-next-btn",
                                       disabled=current_page == total_pages-1, size="sm", color="secondary")
                        ], className="float-end")
                    ], width=6)
                ], className="mb-3")
            ])

            if viz_type == "umap":
                viz_content = create_umap_accessibility_plot(current_scatac_api, current_genes, sample_id, colorscale, threshold, species_key)
            elif viz_type == "cluster":
                viz_content = create_scatac_cluster_umap_plot(current_scatac_api, sample_id, cluster_column)
            elif viz_type == "stats":
                viz_content = create_accessibility_stats_plot(current_scatac_api, current_genes, sample_id, cluster_column, species_key)
            else:
                viz_content = dbc.Alert(f"Unknown visualization type: {viz_type}", color="warning", className="mt-4")

            return html.Div([
                pagination_controls,
                viz_content
            ])

        except Exception as e:
            logger.error(f"❌ Error creating scATAC visualization: {e}")
            traceback.print_exc()
            return dbc.Alert(f"Error creating visualization: {str(e)}", color="danger", className="mt-4")

    # --- scATAC Pagination Callbacks ---
    @app.callback(
        Output("scatac-page-store", "data"),
        [Input("scatac-prev-btn", "n_clicks"),
         Input("scatac-next-btn", "n_clicks")],
        [State("scatac-page-store", "data"),
         State("target-gene-list-store", "data")]
    )
    def update_scatac_pagination(prev_clicks, next_clicks, current_page, target_gene_ids):
        """Handle scATAC pagination controls"""
        ctx = callback_context

        if not target_gene_ids:
            return 0

        current_page = current_page or 0
        genes_per_page = 4
        total_genes = len(target_gene_ids)
        total_pages = (total_genes + genes_per_page - 1) // genes_per_page

        if ctx.triggered:
            button_id = ctx.triggered[0]['prop_id'].split('.')[0]
            if button_id == "scatac-prev-btn" and current_page > 0:
                current_page -= 1
            elif button_id == "scatac-next-btn" and current_page < total_pages - 1:
                current_page += 1

        return current_page

    @app.callback(
        Output("scatac-page-store", "data", allow_duplicate=True),
        [Input("target-gene-list-store", "data")],
        prevent_initial_call=True
    )
    def reset_scatac_page(target_gene_ids):
        """Reset scATAC page to 0 when gene list changes"""
        return 0
