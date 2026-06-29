"""Callbacks for single-cell RNA-seq (scRNA) analysis."""

from __future__ import annotations

import logging
import traceback

import dash
import dash_bootstrap_components as dbc
from dash import Input, Output, State, callback_context, html, no_update

from .apis import get_scrna_api
from .plots import create_cluster_umap_plot, create_umap_plot, create_violin_plot


def register(app, *, logger: logging.Logger) -> None:
    """Register scRNA-seq callbacks."""

    # Study → Sample cascade callback
    @app.callback(
        [Output("scrna-sample-select", "options"),
         Output("scrna-sample-select", "value")],
        Input("scrna-study-select", "value"),
        State("current-species-store", "data")
    )
    def update_scrna_sample_options(study_id, species_key):
        """Update sample options based on selected study"""
        current_scrna_api = get_scrna_api(species_key)
        if not study_id or current_scrna_api is None:
            return [], None

        samples = current_scrna_api.get_samples_by_study(study_id)
        options = []
        for sample in samples:
            n_cells = sample['metadata'].get('n_cells', 0)
            assay_used = sample['metadata'].get('assay_used', 'RNA')

            # Label only shows assay (when not default) and cell count; data version is shown elsewhere
            assay_label = f" [{assay_used}]" if assay_used != 'RNA' else ""
            label = f"{sample['display_name']} ({n_cells:,} cells){assay_label}"
            options.append({'label': label, 'value': sample['id']})

        default_value = samples[0]['id'] if samples else None

        return options, default_value

    @app.callback(
        [Output("scrna-data-version-radio", "options"),
         Output("scrna-data-version-radio", "value")],
        Input("scrna-sample-select", "value"),
        State("scrna-data-version-radio", "value"),
        State("current-species-store", "data")
    )
    def update_scrna_version_controls(sample_id, current_value, species_key):
        """Update data version selector options and selected value."""
        base_options = [
            {'label': '📊 Normalized log1p', 'value': 'norm'},
            {'label': '🪄 MAGIC Smoothed', 'value': 'magic'}
        ]

        current_scrna_api = get_scrna_api(species_key)
        if not sample_id or current_scrna_api is None:
            # Disable options when sample is not ready
            default_options = [
                {**opt, 'disabled': True} for opt in base_options
            ]
            return default_options, 'norm'

        try:
            available_versions = current_scrna_api.get_available_data_versions(sample_id) or {}
            has_magic = 'magic' in available_versions
            has_norm = 'norm' in available_versions
            # Treat legacy datasets (without metadata) as having raw/default data available
            has_raw = 'raw' in available_versions or 'default' in available_versions or not available_versions

            options = []
            for opt in base_options:
                value = opt['value']
                disabled = False
                if value == 'magic':
                    disabled = not has_magic
                elif value == 'norm':
                    disabled = not has_norm
                options.append({**opt, 'disabled': disabled})

            allowed_values = [opt['value'] for opt in options if not opt.get('disabled')]
            if current_value in allowed_values:
                selected_value = current_value
            elif 'norm' in allowed_values:
                selected_value = 'norm'
            elif 'magic' in allowed_values:
                selected_value = 'magic'
            else:
                selected_value = allowed_values[0] if allowed_values else 'norm'

            return options, selected_value

        except Exception as e:
            options = [{**opt, 'disabled': True} for opt in base_options]
            return options, 'norm'

    # Sample → Cluster method cascade callback
    @app.callback(
        [Output("scrna-cluster-select", "options"),
         Output("scrna-cluster-select", "value")],
        Input("scrna-sample-select", "value"),
        State("current-species-store", "data")
    )
    def update_scrna_cluster_options(sample_id, species_key):
        """Update cluster method options based on selected sample"""
        current_scrna_api = get_scrna_api(species_key)
        if not sample_id or current_scrna_api is None:
            return [], None

        # Get cluster info from the selected sample
        cluster_info = current_scrna_api.get_cluster_info(sample_id)

        if cluster_info:
            options = [{'label': col.replace('_', ' ').title(), 'value': col} for col in cluster_info.keys()]

            # Prefer options containing "celltype" (or similar), then any "cluster", else first
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
                {'label': 'Seurat Clusters', 'value': 'seurat_clusters'},
                {'label': 'Resolution 0.2', 'value': 'RNA_snn_res.0.2'}
            ]
            default_value = 'seurat_clusters'

        return options, default_value

    # Toggle scRNA control states based on visualization type
    @app.callback(
        [Output("scrna-cluster-select", "disabled"),
         Output("scrna-colorscale-select", "disabled"),
         Output("scrna-threshold-slider", "disabled")],
        Input("scrna-viz-type-radio", "value"),
        Input("scrna-sample-select", "value")
    )
    def toggle_scrna_control_states(viz_type, sample_id):
        """Enable or disable controls depending on the selected visualization"""
        has_sample = bool(sample_id)
        cluster_disabled = (viz_type not in ("cluster", "violin")) or not has_sample
        color_disabled = (viz_type != "umap") or not has_sample
        threshold_disabled = (viz_type != "umap") or not has_sample
        return cluster_disabled, color_disabled, threshold_disabled

    @app.callback(
        Output("scrna-content", "children"),
        [Input("target-gene-list-store", "data"),
         Input("scrna-sample-select", "value"),
         Input("scrna-cluster-select", "value"),
         Input("scrna-viz-type-radio", "value"),
         Input("scrna-colorscale-select", "value"),
         Input("scrna-data-version-radio", "value"),
         Input("scrna-threshold-slider", "value"),
         Input("scrna-page-store", "data"),
         Input("main-tabs", "value")],
        State("current-species-store", "data")
    )
    def update_scrna_content(target_gene_ids, sample_id, cluster_column, viz_type, colorscale, data_version, threshold, current_page, active_tab, species_key):
        """Update scRNA-seq visualization content"""
        logger.info(f"🧬 update_scrna_content called:")
        logger.info(f"  - active_tab: {active_tab}")
        logger.info(f"  - target_gene_ids: {target_gene_ids}")
        logger.info(f"  - sample_id: {sample_id}")
        logger.info(f"  - viz_type: {viz_type}")

        if active_tab != "tab-scrna":
            return no_update

        current_scrna_api = get_scrna_api(species_key)
        if current_scrna_api is None:
            return dbc.Alert("scRNA API is not available.", color="danger", className="mt-4")

        if not sample_id:
            return dbc.Alert("No sample selected.", color="warning", className="mt-4")

        if not target_gene_ids:
            return dbc.Alert("Select genes from the sidebar to visualize single-cell expression.", color="info", className="mt-4")

        # Calculate pagination
        genes_per_page = 4
        current_page = current_page or 0
        start_idx = current_page * genes_per_page
        end_idx = start_idx + genes_per_page
        current_genes = target_gene_ids[start_idx:end_idx]

        logger.info(f"🧬 scRNA Pagination Debug:")
        logger.info(f"  - Total genes: {len(target_gene_ids)}")
        logger.info(f"  - Current page: {current_page}")
        logger.info(f"  - Page range: {start_idx} to {end_idx}")
        logger.info(f"  - Current genes: {current_genes}")

        if not current_genes:
            return dbc.Alert("No genes to display on this page.", color="warning", className="mt-4")

        try:
            # Create pagination controls (similar to landscape view)
            total_pages = (len(target_gene_ids) + genes_per_page - 1) // genes_per_page

            pagination_controls = html.Div([
                dbc.Row([
                    dbc.Col([
                        html.H6(f"scRNA-seq Visualization", className="mb-0")
                    ], width=6),
                    dbc.Col([
                        dbc.ButtonGroup([
                            dbc.Button("◀ Previous", id="scrna-prev-btn",
                                     disabled=current_page == 0, size="sm", color="secondary"),
                            dbc.Button(f"{current_page + 1} / {total_pages}",
                                     color="light", disabled=True, size="sm"),
                            dbc.Button("Next ▶", id="scrna-next-btn",
                                     disabled=current_page == total_pages-1, size="sm", color="secondary")
                        ], className="float-end")
                    ], width=6)
                ], className="mb-3")
            ])

            # Create visualization content
            if viz_type == "umap":
                viz_content = create_umap_plot(current_scrna_api, current_genes, sample_id, colorscale, threshold, data_version)
            elif viz_type == "cluster":
                viz_content = create_cluster_umap_plot(current_scrna_api, sample_id, cluster_column)
            elif viz_type == "violin":
                viz_content = create_violin_plot(current_scrna_api, current_genes, sample_id, cluster_column, data_version)
            else:
                viz_content = dbc.Alert(f"Unknown visualization type: {viz_type}", color="warning", className="mt-4")

            # Return combined content with pagination controls
            return html.Div([
                pagination_controls,
                viz_content
            ])

        except Exception as e:
            logger.error(f"❌ Error creating scRNA visualization: {e}")
            traceback.print_exc()
            return dbc.Alert(f"Error creating visualization: {str(e)}", color="danger", className="mt-4")

    # --- scRNA Pagination Callbacks ---
    @app.callback(
        Output("scrna-page-store", "data"),
        [Input("scrna-prev-btn", "n_clicks"),
         Input("scrna-next-btn", "n_clicks")],
        [State("scrna-page-store", "data"),
         State("target-gene-list-store", "data")]
    )
    def update_scrna_pagination(prev_clicks, next_clicks, current_page, target_gene_ids):
        """Handle scRNA pagination controls"""
        ctx = callback_context

        if not target_gene_ids:
            return 0

        current_page = current_page or 0
        genes_per_page = 4
        total_genes = len(target_gene_ids)
        total_pages = (total_genes + genes_per_page - 1) // genes_per_page  # Ceiling division

        # Determine which button was clicked
        if ctx.triggered:
            button_id = ctx.triggered[0]['prop_id'].split('.')[0]
            if button_id == "scrna-prev-btn" and current_page > 0:
                current_page -= 1
            elif button_id == "scrna-next-btn" and current_page < total_pages - 1:
                current_page += 1

        return current_page

    # Reset scRNA page when gene list changes
    @app.callback(
        Output("scrna-page-store", "data", allow_duplicate=True),
        [Input("target-gene-list-store", "data")],
        prevent_initial_call=True
    )
    def reset_scrna_page(target_gene_ids):
        """Reset page to 0 when gene list changes"""
        return 0
