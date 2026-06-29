import logging
import time
from typing import Any, Dict, List, Optional

import dash_bootstrap_components as dbc
import numpy as np
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
from dash import dcc, html

from .apis import (
    get_gene_search_api,
    get_bulk_multi_api,
    get_landscape_view_api,
    is_bulkmulti_data_available,
    _normalize_species_key,
)
from .species_data import bulkmulti_data_missing_alert


logger = logging.getLogger("miniodp.dash.plots")

def _build_cluster_violin(plot_data: pd.DataFrame, value_key: str, title: str, cluster_order: Optional[List[str]] = None) -> dcc.Graph:
    """Render a Seurat-like violin plot locked at 0 with width scaling."""
    max_val = float(plot_data[value_key].max()) if not plot_data.empty else 0.0
    max_val = max_val if max_val > 0 and np.isfinite(max_val) else 0.0
    y_upper = max_val * 1.05 if max_val > 0 else 1.0
    span_vals = [0, max_val] if max_val > 0 else [0, 1.0]
    if cluster_order is None:
        cluster_order = sorted(pd.unique(plot_data['Cluster'].astype(str)))

    fig = px.violin(
        plot_data,
        x='Cluster',
        y=value_key,
        color='Cluster',
        title=title,
        points=False,
        box=False,
        category_orders={'Cluster': cluster_order},
    )

    fig.update_traces(
        spanmode="hard",
        span=span_vals,
        scalemode="width",
        meanline_visible=False,
        width=0.9,
        line=dict(width=0.6),
    )
    fig.update_yaxes(range=[0, y_upper], zeroline=True, zerolinewidth=1)
    fig.update_layout(
        title_font_size=16,
        margin=dict(l=20, r=20, t=50, b=20),
        plot_bgcolor='white',
        showlegend=False,
        height=350,
        violinmode="group",
        violingap=0,
    )
    return dcc.Graph(figure=fig, style={'width': '100%'})

def describe_version_label(version_key: str, version_info: Dict[str, str]) -> str:
    """Generate human-friendly label for a data version."""
    type_hint = (version_info.get('type') or '').lower()

    type_mapping = {
        'magic_imputed': 'MAGIC smoothed',
        'raw_counts': 'Raw counts',
        'normalized': 'Normalized log1p',
        'unknown': 'Legacy default'
    }

    version_mapping = {
        'magic': 'MAGIC smoothed',
        'raw': 'Raw counts',
        'default': 'Legacy default'
    }

    if type_hint in type_mapping:
        return type_mapping[type_hint]
    if version_key in version_mapping:
        return version_mapping[version_key]
    if version_info.get('description'):
        return version_info['description']
    return version_key.upper()

def convert_gene_ids_to_names(species_key, gene_ids):
    """Convert gene IDs to gene names using species adapter"""
    species = _normalize_species_key(species_key)
    if not species:
        logger.warning("Gene ID conversion requested without species context")
        return gene_ids  # Return original if no species context
        
    try:
        current_gene_search_api = get_gene_search_api(species)
        if current_gene_search_api:
            adapter = current_gene_search_api.adapter
        else:
            from backend.species_adapters import create_species_adapter
            adapter = create_species_adapter(species)
        
        # Use adapter's convert method
        return adapter.convert_ids_to_gene_names(gene_ids)
        
    except Exception as e:
        logger.error(f"Error converting gene IDs: {e}")
        return gene_ids  # Return original if conversion fails

def create_umap_plot(scrna_api, gene_ids, sample_id, colorscale, threshold, data_version='auto'):
    """Create embedding feature plot for gene expression (auto-detects UMAP/tSNE/PCA)"""
    t0 = time.perf_counter()
    try:
        def _wrap_graph(fig, notice_text=None, notice_color="secondary"):
            if not notice_text:
                return dcc.Graph(figure=fig)
            return html.Div([
                dbc.Alert(notice_text, color=notice_color, className="mb-1", fade=False, dismissable=True),
                dcc.Graph(figure=fig)
            ])
        # Get embedding coordinates with automatic method detection
        t_emb_start = time.perf_counter()
        coord_x, coord_y, method_name = scrna_api.get_embedding_coordinates(sample_id)
        t_emb = time.perf_counter() - t_emb_start
        if len(coord_x) == 0 or len(coord_y) == 0:
            return dbc.Alert(f"Embedding coordinates not available for this sample.", color="warning")
        
        # Get cell metadata for cluster information
        t_obs_start = time.perf_counter()
        obs = scrna_api.get_cell_metadata(sample_id)
        t_obs = time.perf_counter() - t_obs_start
        if obs.empty:
            return dbc.Alert("Cell metadata not available.", color="warning")

        t_version_start = time.perf_counter()
        resolved_version, resolved_info = scrna_api.determine_data_version(sample_id, data_version)
        version_label = describe_version_label(resolved_version, resolved_info) if resolved_info else describe_version_label(resolved_version, {})
        t_version = time.perf_counter() - t_version_start

        # Create plots for each gene
        plots = []
        logger.info(f"🎯 {method_name}: Creating plots for {len(gene_ids)} genes: {gene_ids}")
        t_expr_start = time.perf_counter()
        expr_df, gene_mapping, display_labels = scrna_api.get_gene_expression_batch(gene_ids, sample_id, data_version)
        t_expr = time.perf_counter() - t_expr_start
        if expr_df.empty:
            logger.warning("⚠️ No expression matrix found for requested genes")
            return dbc.Alert("No expression data available for selected genes.", color="warning")

        t_plot_total = 0.0
        coord_x_arr = np.asarray(coord_x)
        coord_y_arr = np.asarray(coord_y)
        for gene_id in gene_ids:  # Process all genes passed to function preserving order
            try:
                resolved_gene_id = gene_mapping.get(gene_id, gene_id)
                if resolved_gene_id not in expr_df.index:
                    logger.warning(f"⚠️ {method_name}: No expression data for gene: {gene_id} - SKIPPING")
                    placeholder_fig = px.scatter(
                        x=coord_x_arr,
                        y=coord_y_arr,
                        color=np.zeros_like(coord_x_arr),
                        color_continuous_scale=[(0, "#f5f5f5"), (1, "#f5f5f5")],
                        labels={'x': f'{method_name} 1', 'y': f'{method_name} 2'},
                        title=f"{gene_id}",
                        width=600,
                        height=500
                    )
                    placeholder_fig.update_layout(coloraxis_showscale=False, plot_bgcolor='white', margin=dict(l=20, r=20, t=50, b=20))
                    placeholder_fig.update_traces(marker=dict(size=3, opacity=0.45))
                    plots.append(_wrap_graph(placeholder_fig, "No data detected", "secondary"))
                    continue

                expression = expr_df.loc[resolved_gene_id].values
                display_name = display_labels.get(gene_id, resolved_gene_id)
                
                # Apply threshold; if everything is filtered out, set all to zero and warn (English hint, not in title)
                applied_threshold = threshold if threshold is not None else 0
                num_above = np.count_nonzero(expression > applied_threshold)
                expression_max = float(np.max(expression)) if expression.size > 0 else 0.0
                color_cmax = expression_max if expression_max > 0 else 1.0  # keep zero at lowest color
                threshold_warning = None
                if expression_max <= 0:
                    logger.info(f"⚠️ {method_name}: No expression detected for {display_name}, colors set to 0")
                    expression_thresholded = np.zeros_like(expression)
                    threshold_note = ""
                    threshold_warning = "No data detected"
                elif num_above == 0:
                    logger.info(
                        f"⚠️ {method_name}: Threshold {applied_threshold} filters all cells for {display_name}, "
                        "all colors set to 0"
                    )
                    expression_thresholded = np.zeros_like(expression)
                    threshold_note = ""
                    threshold_warning = "Current threshold filters all cells."
                else:
                    expression_thresholded = np.where(expression > applied_threshold, expression, 0)
                    threshold_note = ""

                # Draw low-expression cells first to reduce overplotting
                sort_idx = np.argsort(expression_thresholded)
                x_sorted = coord_x_arr[sort_idx]
                y_sorted = coord_y_arr[sort_idx]
                color_sorted = expression_thresholded[sort_idx]

                # Normalize colorscale selection
                selected_scale = colorscale
                show_colorbar = True
                if threshold_warning == "No data detected":
                    selected_scale = [(0, "#f5f5f5"), (1, "#f5f5f5")]
                    color_cmax = 1.0
                    show_colorbar = False
                else:
                    if not selected_scale:
                        selected_scale = px.colors.sequential.Oranges
                    elif isinstance(selected_scale, str):
                        if selected_scale.lower() == 'ylgnbu':
                            selected_scale = px.colors.sequential.YlGnBu
                        elif selected_scale.lower() == 'viridis':
                            selected_scale = px.colors.sequential.Viridis
                        elif selected_scale.lower() == 'oranges':
                            selected_scale = px.colors.sequential.Oranges

                # Create scatter plot
                t_fig_start = time.perf_counter()
                fig = px.scatter(
                    x=x_sorted, 
                    y=y_sorted,
                    color=color_sorted,
                    color_continuous_scale=selected_scale,
                    title=f"{display_name} Expression – {method_name}, {version_label}{threshold_note}",
                    labels={'color': 'Expression', 'x': f'{method_name} 1', 'y': f'{method_name} 2'},
                    width=600,
                    height=500
                )
                # Keep zero at lowest color; clamp upper bound to avoid extreme outliers flattening the scale
                safe_cmax = color_cmax
                try:
                    pct99 = float(np.percentile(expression_thresholded, 99)) if expression_thresholded.size > 0 else 0.0
                    if pct99 > 0:
                        safe_cmax = min(color_cmax, pct99)
                except Exception:
                    pass
                if safe_cmax <= 0:
                    safe_cmax = color_cmax if color_cmax > 0 else 1.0

                fig.update_layout(
                    coloraxis_cmin=0,
                    coloraxis_cmax=safe_cmax,
                    coloraxis_showscale=show_colorbar
                )
                fig.update_traces(marker=dict(size=3, opacity=0.7))
                fig.update_layout(
                    title_font_size=16,
                    margin=dict(l=20, r=20, t=50, b=20),
                    plot_bgcolor='white'
                )

                alert_color = "secondary" if threshold_warning == "No data detected" else "warning"
                graph_component = _wrap_graph(fig, threshold_warning, alert_color)

                plots.append(graph_component)
                t_plot_total += time.perf_counter() - t_fig_start
                
            except Exception as e:
                logger.warning(f"⚠️ Error processing gene {gene_id}: {e}")
                continue
        
        if not plots:
            return dbc.Alert("No expression data found for selected genes.", color="warning")
        
        logger.info(
            "⏱️ UMAP timing dataset=%s genes=%s version=%s | emb=%.2fs obs=%.2fs select_ver=%.2fs expr=%.2fs plot=%.2fs total=%.2fs",
            sample_id,
            gene_ids,
            resolved_version,
            t_emb,
            t_obs,
            t_version,
            t_expr,
            t_plot_total,
            time.perf_counter() - t0,
        )
        logger.info(f"✅ {method_name}: Successfully created {len(plots)} plots out of {len(gene_ids)} genes")
        
        # Arrange plots in grid
        if len(plots) == 1:
            return dbc.Row([dbc.Col(plots[0], width=6)])
        elif len(plots) == 2:
            return dbc.Row([
                dbc.Col(plots[0], width=6),
                dbc.Col(plots[1], width=6)
            ])
        else:
            return html.Div([
                dbc.Row([
                    dbc.Col(plots[0], width=6),
                    dbc.Col(plots[1], width=6)
                ]),
                dbc.Row([
                    dbc.Col(plots[i], width=6) for i in range(2, min(len(plots), 4))
                ]) if len(plots) > 2 else None
            ])
        
    except Exception as e:
        logger.error(f"❌ Error creating embedding plot: {e}")
        return dbc.Alert(f"Error creating embedding plot: {str(e)}", color="danger")

def create_cluster_umap_plot(scrna_api, sample_id, cluster_column):
    """Create UMAP plot colored by cluster labels"""
    try:
        coord_x, coord_y, method_name = scrna_api.get_embedding_coordinates(sample_id)
        if len(coord_x) == 0 or len(coord_y) == 0:
            return dbc.Alert("Embedding coordinates not available for this sample.", color="warning")

        obs = scrna_api.get_cell_metadata(sample_id)
        if obs.empty or cluster_column not in obs.columns:
            return dbc.Alert(f"Cluster column '{cluster_column}' not found.", color="warning")

        clusters = obs[cluster_column].astype(str)

        fig = px.scatter(
            x=coord_x,
            y=coord_y,
            color=clusters,
            title=f"Cell Clusters ({method_name})",
            labels={'color': 'Cluster', 'x': f'{method_name} 1', 'y': f'{method_name} 2'},
            width=800,
            height=600
        )

        fig.update_traces(marker=dict(size=4, opacity=0.7))
        fig.update_layout(
            title_font_size=16,
            margin=dict(l=20, r=20, t=50, b=20),
            plot_bgcolor='white',
            legend=dict(title="Cluster", itemsizing='constant')
        )

        return dcc.Graph(figure=fig)

    except Exception as e:
        logger.error(f"❌ Error creating cluster UMAP plot: {e}")
        return dbc.Alert(f"Error creating cluster plot: {str(e)}", color="danger")

def create_violin_plot(scrna_api, gene_ids, sample_id, cluster_column, data_version='auto'):
    """Create violin plot for gene expression by clusters"""
    try:
        # Get raw expression data and cell metadata
        expr_df, gene_mapping, display_labels = scrna_api.get_gene_expression_batch(gene_ids, sample_id, data_version)
        obs = scrna_api.get_cell_metadata(sample_id)
        
        if expr_df.empty or obs.empty or cluster_column not in obs.columns:
            return dbc.Alert("No expression or cluster data available.", color="warning")
        
        # Align cell names
        common_cells = expr_df.columns.intersection(obs.index)
        if len(common_cells) == 0:
            return dbc.Alert("No overlapping cells between expression and metadata.", color="warning")
        expr_subset = expr_df[common_cells]
        obs_subset = obs.loc[common_cells]
        clusters = obs_subset[cluster_column].astype(str)
        
        # Create violin plots
        plots = []
        logger.info(f"🎯 VIOLIN: Creating plots for {len(gene_ids)} genes: {gene_ids}")

        resolved_version, resolved_info = scrna_api.determine_data_version(sample_id, data_version)
        version_label = describe_version_label(resolved_version, resolved_info) if resolved_info else describe_version_label(resolved_version, {})

        plotted_genes = set()
        for gene_id in gene_ids:  # Process all genes passed to function
            try:
                resolved_gene_id = gene_mapping.get(gene_id, gene_id)
                if resolved_gene_id in plotted_genes:
                    logger.info(f"⚠️ VIOLIN: Duplicate resolved gene {resolved_gene_id}, skipping duplicate plot")
                    continue
                if resolved_gene_id not in expr_subset.index:
                    logger.warning(f"⚠️ VIOLIN: No expression data for gene: {gene_id} - SKIPPING")
                    continue

                if expr_subset.empty or obs_subset.empty:
                    continue

                expression_values = expr_subset.loc[resolved_gene_id]
                plotted_genes.add(resolved_gene_id)
                
                # Create DataFrame for plotting
                plot_data = pd.DataFrame({
                    'Expression': expression_values,
                    'Cluster': clusters
                })
                
                display_name = display_labels.get(gene_id, resolved_gene_id)
                title = f"{display_name} Expression by Cluster ({version_label})"
                plots.append(_build_cluster_violin(plot_data, 'Expression', title))
                
            except Exception as e:
                logger.warning(f"⚠️ Error processing gene {gene_id}: {e}")
                continue
        
        if not plots:
            return dbc.Alert("No expression data found for selected genes.", color="warning")
        
        logger.info(f"✅ VIOLIN: Successfully created {len(plots)} plots out of {len(gene_ids)} genes")
        
        # Single column layout - each gene takes full width
        rows = [dbc.Row([dbc.Col(plot, width=12)], className="mb-2") for plot in plots]
        return html.Div(rows)
        
    except Exception as e:
        logger.error(f"❌ Error creating violin plot: {e}")
        return dbc.Alert(f"Error creating violin plot: {str(e)}", color="danger")

def create_umap_accessibility_plot(scatac_api, gene_ids, dataset_id, colorscale, threshold, species_key):
    """Create UMAP plot for gene activity (GAS)"""
    try:
        def _wrap_graph(fig, notice_text=None, notice_color="secondary"):
            if not notice_text:
                return dcc.Graph(figure=fig)
            return html.Div([
                dbc.Alert(notice_text, color=notice_color, className="mb-1", fade=False, dismissable=True),
                dcc.Graph(figure=fig)
            ])
        umap_x, umap_y = scatac_api.get_umap_coordinates(dataset_id)
        if len(umap_x) == 0 or len(umap_y) == 0:
            return dbc.Alert("UMAP coordinates not available for this dataset.", color="warning")

        obs = scatac_api.get_cell_metadata(dataset_id)
        if obs.empty:
            return dbc.Alert("Cell metadata not available.", color="warning")
        activity_df, gene_mapping, display_labels = scatac_api.get_gene_activity_batch(gene_ids, dataset_id)
        if activity_df.empty:
            return dbc.Alert("No gene activity available for selected genes.", color="warning")

        obs_index = pd.Index(obs.index, dtype=str)
        if not activity_df.columns.equals(obs_index):
            # Require overlap, but prefer obs order
            common = obs_index.intersection(activity_df.columns)
            if common.empty:
                return dbc.Alert("Gene activity matrix and metadata share no overlapping cells.", color="warning")
            activity_df = activity_df.reindex(columns=obs_index)

        plotted_genes = set()
        plots = []
        method_name = "UMAP"
        umap_x_arr = np.asarray(umap_x)
        umap_y_arr = np.asarray(umap_y)
        for gene_id in gene_ids:
            try:
                resolved_gene_id = gene_mapping.get(gene_id, gene_id)
                if resolved_gene_id not in activity_df.index:
                    placeholder_fig = px.scatter(
                        x=umap_x_arr,
                        y=umap_y_arr,
                        color=np.zeros_like(umap_x_arr),
                        color_continuous_scale=[(0, "#f5f5f5"), (1, "#f5f5f5")],
                        labels={'x': f'{method_name} 1', 'y': f'{method_name} 2'},
                        title=f"{gene_id}",
                        width=600,
                        height=500
                    )
                    placeholder_fig.update_layout(coloraxis_showscale=False, plot_bgcolor='white', margin=dict(l=20, r=20, t=50, b=20))
                    placeholder_fig.update_traces(marker=dict(size=3, opacity=0.45))
                    plots.append(_wrap_graph(placeholder_fig, "No data detected", "secondary"))
                    continue

                activity = np.nan_to_num(activity_df.loc[resolved_gene_id].values, nan=0.0)
                display_name = display_labels.get(gene_id, resolved_gene_id)
                plotted_genes.add(resolved_gene_id)

                applied_threshold = threshold if threshold is not None else 0
                num_above = np.count_nonzero(activity > applied_threshold)
                activity_max = float(np.max(activity)) if activity.size > 0 else 0.0
                color_cmax = activity_max if activity_max > 0 else 1.0
                threshold_warning = None
                if activity_max <= 0:
                    activity_thresholded = np.zeros_like(activity)
                    threshold_warning = "No data detected"
                elif num_above == 0:
                    activity_thresholded = np.zeros_like(activity)
                    threshold_warning = "Current threshold filters all cells."
                else:
                    activity_thresholded = np.where(activity > applied_threshold, activity, 0)

                sort_idx = np.argsort(activity_thresholded)
                x_sorted = umap_x_arr[sort_idx]
                y_sorted = umap_y_arr[sort_idx]
                color_sorted = activity_thresholded[sort_idx]

                selected_scale = colorscale
                show_colorbar = True
                if threshold_warning == "No data detected":
                    selected_scale = [(0, "#f5f5f5"), (1, "#f5f5f5")]
                    color_cmax = 1.0
                    show_colorbar = False
                else:
                    if not selected_scale:
                        selected_scale = px.colors.sequential.Oranges
                    elif isinstance(selected_scale, str):
                        if selected_scale.lower() == 'ylgnbu':
                            selected_scale = px.colors.sequential.YlGnBu
                        elif selected_scale.lower() == 'viridis':
                            selected_scale = px.colors.sequential.Viridis
                        elif selected_scale.lower() == 'oranges':
                            selected_scale = px.colors.sequential.Oranges

                fig = px.scatter(
                    x=x_sorted, 
                    y=y_sorted,
                    color=color_sorted,
                    color_continuous_scale=selected_scale,
                    title=f"{display_name} Gene activity ({method_name}, normalized log1p)",
                    labels={'color': 'Gene activity', 'x': f'{method_name} 1', 'y': f'{method_name} 2'},
                    width=600,
                    height=500
                )

                # Clamp color axis to reduce outlier flattening
                safe_cmax = color_cmax
                try:
                    pct99 = float(np.percentile(activity_thresholded, 99)) if activity_thresholded.size > 0 else 0.0
                    if pct99 > 0:
                        safe_cmax = min(color_cmax, pct99)
                except Exception:
                    pass
                if safe_cmax <= 0:
                    safe_cmax = color_cmax if color_cmax > 0 else 1.0

                fig.update_layout(
                    coloraxis_cmin=0,
                    coloraxis_cmax=safe_cmax,
                    coloraxis_showscale=show_colorbar
                )
                fig.update_traces(marker=dict(size=3, opacity=0.7))
                fig.update_layout(
                    title_font_size=16,
                    margin=dict(l=20, r=20, t=50, b=20),
                    plot_bgcolor='white'
                )

                alert_color = "secondary" if threshold_warning == "No data detected" else "warning"
                graph_component = _wrap_graph(fig, threshold_warning, alert_color)

                plots.append(graph_component)

            except Exception as e:
                logger.warning(f"⚠️ Error processing gene {gene_id}: {e}")
                continue

        if not plots:
            return dbc.Alert("No accessibility data found for selected genes.", color="warning")

        rows = []
        for i in range(0, len(plots), 2):
            row_plots = plots[i:i+2]
            rows.append(dbc.Row([dbc.Col(p, width=6) for p in row_plots], className="mb-4"))
        return html.Div(rows)
        
    except Exception as e:
        logger.error(f"❌ Error creating UMAP accessibility plot: {e}")
        return dbc.Alert(f"Error creating UMAP plot: {str(e)}", color="danger")

def create_scatac_cluster_umap_plot(scatac_api, dataset_id, cluster_column):
    """Create cluster-level UMAP for scATAC."""
    try:
        umap_x, umap_y = scatac_api.get_umap_coordinates(dataset_id)
        if len(umap_x) == 0 or len(umap_y) == 0:
            return dbc.Alert("UMAP coordinates not available for this dataset.", color="warning")

        obs = scatac_api.get_cell_metadata(dataset_id)
        if obs.empty:
            return dbc.Alert("Cell metadata not available.", color="warning")
        if not cluster_column:
            return dbc.Alert("No cluster column selected.", color="warning")
        if cluster_column not in obs.columns:
            return dbc.Alert(f"Cluster column '{cluster_column}' not found.", color="warning")

        clusters = obs[cluster_column].astype(str)
        fig = px.scatter(
            x=umap_x,
            y=umap_y,
            color=clusters,
            title="Cell Clusters (UMAP)",
            labels={'color': 'Cluster', 'x': 'UMAP1', 'y': 'UMAP2'},
            width=800,
            height=600
        )
        fig.update_traces(marker=dict(size=4, opacity=0.7))
        fig.update_layout(
            title_font_size=16,
            margin=dict(l=20, r=20, t=50, b=20),
            plot_bgcolor='white',
            legend=dict(title="Cluster", itemsizing='constant')
        )
        return dcc.Graph(figure=fig)
    except Exception as e:
        logger.error(f"❌ Error creating scATAC cluster UMAP plot: {e}")
        return dbc.Alert(f"Error creating cluster plot: {str(e)}", color="danger")

def create_accessibility_stats_plot(scatac_api, gene_ids, dataset_id, cluster_column, species_key):
    """Create gene activity violin plot by clusters"""
    try:
        obs = scatac_api.get_cell_metadata(dataset_id)
        if obs.empty:
            return dbc.Alert("Cell metadata not available.", color="warning")
        if cluster_column not in obs.columns:
            return dbc.Alert(f"Cluster column {cluster_column} not found.", color="warning")
        clusters = obs[cluster_column].astype(str)

        activity_df, gene_mapping, display_labels = scatac_api.get_gene_activity_batch(gene_ids, dataset_id)
        if activity_df.empty:
            return dbc.Alert("No gene activity available for selected genes.", color="warning")

        # Align cells
        obs_index = pd.Index(obs.index, dtype=str)
        common_cells = obs_index.intersection(activity_df.columns)
        if common_cells.empty:
            return dbc.Alert("Gene activity matrix and metadata share no overlapping cells.", color="warning")
        
        activity_df = activity_df.reindex(columns=common_cells)
        clusters_aligned = clusters.reindex(index=common_cells)

        # Create individual plots for each gene (same style as scRNA)
        plots = []
        plotted_genes = set()
        
        for gene_id in gene_ids:
            resolved_gene_id = gene_mapping.get(gene_id, gene_id)
            if resolved_gene_id in plotted_genes:
                continue
            if resolved_gene_id not in activity_df.index:
                continue
            
            plotted_genes.add(resolved_gene_id)
            display_name = display_labels.get(gene_id, resolved_gene_id)
            
            # Get activity values for this gene
            activity_values = activity_df.loc[resolved_gene_id]
            
            # Create DataFrame for plotting
            plot_data = pd.DataFrame({
                'Activity': activity_values.values,
                'Cluster': clusters_aligned.values
            })
            
            title = f"{display_name} Gene Activity by Cluster"
            plots.append(_build_cluster_violin(plot_data, 'Activity', title))
            
            # Limit to 4 genes max
            if len(plots) >= 4:
                break
        
        if not plots:
            return dbc.Alert("No gene activity data found for selected genes.", color="warning")
        
        # Single column layout - each gene takes full width
        rows = [dbc.Row([dbc.Col(plot, width=12)], className="mb-2") for plot in plots]
        return html.Div(rows)
        
    except Exception as e:
        logger.error(f"❌ Error creating accessibility stats plot: {e}")
        return dbc.Alert(f"Error creating statistics plot: {str(e)}", color="danger")

def _build_bulkmulti_view(species_key, target_gene_ids, view_type, selected_dataset, value_threshold):
    """Render BulkMulti content based on requested options."""
    bulk_multi_api = get_bulk_multi_api(species_key, selected_dataset)
    landscape_view_api = get_landscape_view_api(species_key, selected_dataset)
    
    if not bulk_multi_api:
        return dbc.Alert("BulkMulti API is not initialized for this species.", color="warning")
    if not is_bulkmulti_data_available(species_key, selected_dataset):
        return bulkmulti_data_missing_alert()

    try:
        gene_names = convert_gene_ids_to_names(species_key, target_gene_ids)
        logger.info(f"🧬 Converted gene names: {gene_names}")

        if not gene_names:
            return dbc.Alert("No valid gene names found for BulkMulti analysis", color="warning")

        if view_type == "statistics":
            logger.info("🔄 Creating statistics view (deferred render)...")
            stats = bulk_multi_api.get_linkage_statistics(gene_names)
            return create_bulkmulti_statistics_view(bulk_multi_api, gene_names, value_threshold, stats)

        logger.info("🔄 Creating landscape view (deferred render)...")
        landscape_config = {
            'upstream_kb': 50,
            'downstream_kb': 50,
            'show_els_guide_lines': False,
            'atac_height_ratio': 0.2,
            'value_threshold': value_threshold
        }
        logger.info(f"🎛️ Using default landscape config: {landscape_config}")
        return create_bulkmulti_landscape_view(
            bulk_multi_api,
            landscape_view_api,
            gene_names,
            landscape_config['value_threshold'],
            species_key,
            landscape_config
        )

    except Exception as exc:
        logger.error(f"❌ Error building BulkMulti view: {exc}")
        import traceback
        traceback.print_exc()
        return dbc.Alert(f"Error creating BulkMulti visualization: {str(exc)}", color="danger")

def create_bulkmulti_statistics_view(bulk_api, gene_names, value_threshold, stats, chart_types=None):
    """Create BulkMulti statistics view"""
    if chart_types is None:
        chart_types = ['bar_chart', 'scatter_plot', 'stats_table']
    try:
        if not stats:
            return dbc.Alert("No linkage data found for selected genes", color="info")

        notices = []
        missing_genes = [g for g in gene_names if g not in stats or stats[g].get('total_linkages', 0) == 0]
        if missing_genes:
            notices.append(
                dbc.Alert(
                    f"No linkage data for: {', '.join(missing_genes)} (threshold {value_threshold})",
                    color="warning",
                    className="mb-2",
                    fade=False,
                    dismissable=False
                )
            )
        
        # Create statistics table
        stats_data = []
        for gene_name, gene_stats in stats.items():
            stats_data.append({
                'Gene': gene_name,
                'Total Linkages': gene_stats['total_linkages'],
                'Strong (≥0.5)': gene_stats['strong_linkages'],
                'Moderate (0.2-0.5)': gene_stats['moderate_linkages'],
                'Weak (<0.2)': gene_stats['weak_linkages'],
                'Max Strength': f"{gene_stats['max_linkage_value']:.3f}",
                'Mean Strength': f"{gene_stats['mean_linkage_value']:.3f}",
                'Median Distance': f"{gene_stats['median_distance']/1000:.1f}kb",
                'Chromosomes': ', '.join(gene_stats['chromosomes'][:3])  # Show first 3 chromosomes
            })
        
        stats_df = pd.DataFrame(stats_data)
        
        # Create bar chart for linkage counts
        fig_counts = px.bar(
            stats_df, 
            x='Gene', 
            y=['Strong (≥0.5)', 'Moderate (0.2-0.5)', 'Weak (<0.2)'],
            title="Peak2Gene Linkage Counts by Strength",
            labels={'value': 'Number of Linkages', 'variable': 'Linkage Strength'},
            height=400
        )
        fig_counts.update_layout(
            xaxis_title="Gene",
            yaxis_title="Number of Linkages",
            legend_title="Linkage Strength"
        )
        
        # Create scatter plot for linkage strength vs distance
        linkages_by_gene = bulk_api.get_gene_linkages_by_gene(gene_names, value_threshold)
        linkage_data = []
        for gene_name in gene_names:
            gene_linkages = linkages_by_gene.get(gene_name)
            if gene_linkages is None or gene_linkages.empty:
                continue
            for _, row in gene_linkages.iterrows():
                linkage_data.append({
                    'Gene': gene_name,
                    'Distance (kb)': row['Dist'] / 1000,
                    'Linkage Strength': row['Value'],
                    'Chromosome': row['Chromosome']
                })
        
        fig_scatter = None
        if linkage_data:
            linkage_df = pd.DataFrame(linkage_data)
            fig_scatter = px.scatter(
                linkage_df,
                x='Distance (kb)',
                y='Linkage Strength',
                color='Gene',
                title="Peak2Gene Linkage Strength vs Distance",
                height=400,
                hover_data=['Chromosome']
            )
            fig_scatter.update_layout(
                xaxis_title="Distance from Gene (kb)",
                yaxis_title="Linkage Strength"
            )
        else:
            notices.append(
                dbc.Alert(
                    f"No linkage scatter data at threshold {value_threshold}",
                    color="info",
                    className="mb-2",
                    fade=False,
                    dismissable=False
                )
            )
        
        # Build layout based on selected chart types
        components = notices.copy()
        
        # Statistics table
        if 'stats_table' in chart_types:
            components.append(
                dbc.Row([
                    dbc.Col([
                        html.H5("Linkage Statistics", className="mb-3"),
                        dbc.Table.from_dataframe(stats_df, striped=True, bordered=True, hover=True, size="sm")
                    ], width=12)
                ], className="mb-4")
            )
        
        # Charts row
        chart_cols = []
        if 'bar_chart' in chart_types:
            chart_cols.append(
                dbc.Col([
                    dcc.Graph(figure=fig_counts)
                ], width=6 if 'scatter_plot' in chart_types else 12)
            )
        
        if 'scatter_plot' in chart_types and fig_scatter is not None:
            chart_cols.append(
                dbc.Col([
                    dcc.Graph(figure=fig_scatter)
                ], width=6 if 'bar_chart' in chart_types else 12)
            )
        
        if chart_cols:
            components.append(dbc.Row(chart_cols))
        
        if not components:
            components.append(
                dbc.Alert("No chart types selected. Please choose at least one visualization option.", color="info")
            )
        
        return html.Div(components)
        
    except Exception as e:
        logger.error(f"❌ Error creating BulkMulti statistics view: {e}")
        return dbc.Alert(f"Error creating statistics view: {str(e)}", color="danger")

def create_single_gene_landscape_view(bulk_api, landscape_api, gene_name, page_index, total_genes, landscape_config, value_threshold):
    """Create landscape view for a single gene with pagination info"""
    try:
        notice = None
        placeholder_fig = go.Figure()
        placeholder_fig.update_layout(
            height=600,
            plot_bgcolor="white",
            paper_bgcolor="white",
            xaxis=dict(visible=False),
            yaxis=dict(visible=False),
            margin=dict(l=20, r=20, t=40, b=20)
        )

        # Get single gene information
        gene_info = bulk_api.get_gene_info([gene_name])
        if gene_info.empty:
            notice = dbc.Alert(f"Gene {gene_name} not found", color="warning", className="mb-2", fade=False, dismissable=False)
            return create_paginated_landscape_layout(
                gene_name, page_index, total_genes, placeholder_fig, landscape_config, 
                "NA", 0, 0, notice=notice
            )
        
        gene = gene_info.iloc[0]
        chromosome = gene['Chromosome']
        
        # Calculate region with user-defined extensions
        upstream_bp = landscape_config['upstream_kb'] * 1000
        downstream_bp = landscape_config['downstream_kb'] * 1000
        region_start = gene['Start'] - upstream_bp
        region_end = gene['End'] + downstream_bp
        
        # Create landscape view API configuration
        custom_config = {
            'atac_height_ratio': landscape_config['atac_height_ratio'],
            'show_els_guide_lines': landscape_config['show_els_guide_lines'],
            'sample_order': []  # Use default order
        }
        
        if not landscape_api:
            raise ValueError("LandscapeViewAPI is not initialized")

        # Create the landscape view plot for single gene
        landscape_fig = landscape_api.create_landscape_view_plot(
            gene_names=[gene_name],  # Only current gene
            chromosome=chromosome,
            region_start=region_start,
            region_end=region_end,
            bulk_multi_api=bulk_api,
            value_threshold=value_threshold,
            render_config=custom_config
        )
        if landscape_fig is None or (hasattr(landscape_fig, "data") and len(landscape_fig.data) == 0):
            notice = dbc.Alert(
                f"No landscape data for {gene_name} at threshold {value_threshold}",
                color="info",
                className="mb-2",
                fade=False,
                dismissable=False
            )
            landscape_fig = placeholder_fig
        
        # Create paginated layout
        return create_paginated_landscape_layout(
            gene_name, page_index, total_genes, landscape_fig, landscape_config, 
            chromosome, region_start, region_end, notice=notice
        )
        
    except FileNotFoundError as missing_err:
        logger.error(f"❌ BulkMulti dataset missing while rendering landscape for {gene_name}: {missing_err}")
        return bulkmulti_data_missing_alert()
    except Exception as e:
        logger.error(f"❌ Error creating single gene landscape for {gene_name}: {e}")
        import traceback
        traceback.print_exc()
        return create_paginated_landscape_layout(
            gene_name, page_index, total_genes, placeholder_fig, landscape_config,
            "NA", 0, 0,
            notice=dbc.Alert(f"Error creating landscape for {gene_name}: {str(e)}", color="danger", className="mb-2", fade=False, dismissable=False)
        )

def create_paginated_landscape_layout(gene_name, page_index, total_genes, landscape_fig, 
                                    landscape_config, chromosome, region_start, region_end, notice=None):
    """Create the paginated landscape layout with navigation controls"""
    return html.Div([
        # Pagination controls
        dbc.Card([
            dbc.CardBody([
                dbc.Row([
                    dbc.Col([
                        html.H5(f"Landscape View - Gene {page_index + 1} of {total_genes}", className="mb-1"),
                        html.P(f"Current Gene: {gene_name}", className="text-primary fw-bold mb-0")
                    ], width=6),
                    dbc.Col([
                        dbc.ButtonGroup([
                            dbc.Button("◀ Previous", id="landscape-prev-btn", 
                                     disabled=page_index == 0, size="sm", color="secondary"),
                            dbc.Button(f"{page_index + 1} / {total_genes}", 
                                     color="light", disabled=True, size="sm"),
                            dbc.Button("Next ▶", id="landscape-next-btn", 
                                     disabled=page_index == total_genes-1, size="sm", color="secondary")
                        ], className="float-end")
                    ], width=6)
                ])
            ])
        ], className="mb-3"),
        
        # Gene region info
        dbc.Card([
            dbc.CardBody([
                html.P(f"Showing gene on {chromosome}", className="text-muted mb-1"),
                html.P(f"Region: {int(region_start):,} - {int(region_end):,} bp", className="text-muted small mb-1"),
                html.P(f"Extension: -{landscape_config['upstream_kb']}kb / +{landscape_config['downstream_kb']}kb | ELS Lines: {'On' if landscape_config['show_els_guide_lines'] else 'Off'}", className="text-muted small mb-0")
            ])
        ], className="mb-4"),
        
        # Main landscape view plot
        dbc.Card([
            dbc.CardBody([
                notice if notice is not None else None,
                dcc.Graph(
                    id="landscape-view-plot",
                    figure=landscape_fig,
                    style={'height': '1200px'}
                )
            ])
        ]),
        
        # Note: Store components are managed at the top level to avoid ID conflicts
    ])

def create_bulkmulti_landscape_view(bulk_api, landscape_api, gene_names, value_threshold, species_key, landscape_config=None):
    """Create BulkMulti landscape view with pagination support"""
    if landscape_config is None:
        landscape_config = {
            'upstream_kb': 50,
            'downstream_kb': 50, 
            'show_els_guide_lines': False,
            'max_linkages_display': 15,
            'atac_height_ratio': 0.2
        }
    try:
        if not gene_names or not bulk_api or not landscape_api:
            return dbc.Alert("No genes selected or required APIs not available", color="warning")
        if not is_bulkmulti_data_available(species_key, getattr(bulk_api, "dataset", None)):
            return bulkmulti_data_missing_alert()
        
        # Add value_threshold to landscape_config for consistency
        landscape_config['value_threshold'] = value_threshold
        
        # If only one gene, display directly without pagination
        if len(gene_names) == 1:
            return create_single_gene_landscape_view(bulk_api, landscape_api, gene_names[0], 0, 1, landscape_config, value_threshold)
        
        # Multiple genes: display first gene with pagination
        return create_paginated_landscape_view_initial(bulk_api, landscape_api, gene_names, landscape_config, value_threshold)
        
    except Exception as e:
        logger.error(f"❌ Error creating BulkMulti landscape view: {e}")
        return dbc.Alert(f"Error creating landscape view: {str(e)}", color="danger")

def create_paginated_landscape_view_initial(bulk_api, landscape_api, gene_names, landscape_config, value_threshold):
    """Create initial paginated landscape view showing first gene"""
    try:
        # Display first gene with pagination controls
        return html.Div([
            # Store gene list and pagination state
            dcc.Store(id="landscape-genes-store", data=gene_names),
            dcc.Store(id="landscape-config-store", data=landscape_config),
            dcc.Store(id="landscape-threshold-store", data=value_threshold),
            dcc.Store(id="landscape-page-store", data=0),  # Start at page 0
            
            # Create single gene view for first gene
            create_single_gene_landscape_view(bulk_api, landscape_api, gene_names[0], 0, len(gene_names), landscape_config, value_threshold)
        ])
        
    except Exception as e:
        logger.error(f"❌ Error creating initial paginated landscape view: {e}")
        return dbc.Alert(f"Error creating paginated landscape view: {str(e)}", color="danger")
