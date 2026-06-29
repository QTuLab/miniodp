"""Callbacks for Bulk RNA expression analysis."""

from __future__ import annotations

import logging
from datetime import datetime

import dash
import dash_bootstrap_components as dbc
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
from dash import Input, Output, State, dcc, html, no_update

from .apis import get_expression_api, get_gene_search_api


def register(app, *, logger: logging.Logger) -> None:
    """Register bulk RNA expression callbacks."""

    @app.callback(
        [Output("bulk-rna-content", "children"),
         Output("expression-matrix-store", "data")],
        [Input("target-gene-list-store", "data"),
         Input("expression-sources-checklist", "value"),
         Input("expression-viz-type-radio", "value"),
         Input("expression-transform-radio", "value"),
         Input("expression-cluster-switch", "value"),
         Input("expression-show-values-switch", "value"),
         Input("custom-sample-selection-switch", "value"),
         Input({"type": "sample-checklist", "source": dash.ALL}, "value"),
         Input("main-tabs", "value")],
        State("current-species-store", "data")
    )
    def update_bulk_rna_content(target_gene_ids, selected_sources, viz_type, transform_method,
                               cluster_genes, show_values, custom_sample_selection, selected_samples_by_source, active_tab, species_key):
        logger.info(f"🧬 update_bulk_rna_content called:")
        logger.info(f"  - active_tab: {active_tab}")
        logger.info(f"  - target_gene_ids: {target_gene_ids}")
        logger.info(f"  - selected_sources: {selected_sources}")
        logger.info(f"  - viz_type: {viz_type}")
        logger.info(f"  - transform_method: {transform_method}")
        logger.info(f"  - custom_sample_selection: {custom_sample_selection}")
        logger.info(f"  - species_key: {species_key}")

        if active_tab != "tab-bulk-rna":
            logger.error("  ❌ Not BulkRNA tab, returning no_update")
            return no_update, {}

        if not target_gene_ids:
            return dbc.Alert("Select genes from the sidebar to see their expression patterns.", color="info", className="mt-4"), {}

        # Fetch species-specific expression API
        current_expression_api = get_expression_api(species_key)
        if current_expression_api is None:
            logger.error("  ❌ Expression API is unavailable for species_key=%s", species_key)
            return dbc.Alert("Expression API is not available.", color="danger", className="mt-4"), {}

        if not selected_sources:
            logger.warning("  ⚠️ No sources selected for species_key=%s", species_key)
            return dbc.Alert("Please select at least one study.", color="warning", className="mt-4"), {}

        try:
            # Determine which samples to use
            selected_samples = None
            if custom_sample_selection and selected_samples_by_source:
                selected_samples = []
                for samples_list in selected_samples_by_source:
                    if not samples_list:
                        continue
                    selected_samples.extend([
                        sample_key for sample_key in samples_list if isinstance(sample_key, str)
                    ])
                logger.info(f"  📋 Custom sample selection: {len(selected_samples)} samples")

            # Convert primary_ids to common_ids for expression data queries
            # For medaka: IGDB::ENS -> IGDB; For others: ENS unchanged
            common_gene_ids = current_expression_api.adapter.get_common_ids(target_gene_ids)

            # Get expression matrix
            logger.info(f"  ✅ Querying expression data for {len(common_gene_ids)} common_ids from {len(selected_sources)} sources")
            logger.info(f"  📝 Common gene IDs: {common_gene_ids}")
            logger.info(f"  📝 Current species: {species_key}")
            logger.info(f"  📝 Expression API gene ID column: {current_expression_api.gene_id_column}")
            matrix, metadata = current_expression_api.get_expression_matrix(
                gene_ids=common_gene_ids,
                sources=selected_sources,
                samples=selected_samples
            )

            if matrix.empty:
                return dbc.Alert(f"No expression data found for the selected genes and sources. {metadata.get('error', '')}",
                               color="warning", className="mt-4"), {}

            logger.info(f"  📊 Retrieved matrix: {matrix.shape}")

            # Transform data
            if transform_method != 'raw':
                matrix = current_expression_api.transform_expression_data(matrix, method=transform_method)
                logger.info(f"  🔄 Applied {transform_method} transformation")

            # Expand matrix to show one row per primary_id (for medaka one-to-many IGDB->ENS)
            # Build mapping: primary_id -> common_id (IGDB)
            current_gene_search_api = get_gene_search_api(species_key)
            expanded_rows = []
            expanded_labels = []
            for pid in target_gene_ids:
                parsed = current_expression_api.adapter.parse_primary_id(pid)
                common_id = parsed.get('igdb') or parsed.get('ens') or pid
                if common_id in matrix.index:
                    expanded_rows.append(matrix.loc[common_id])
                    # Get gene name for display
                    gene_name = current_gene_search_api.get_gene_names_by_ids([common_id]).get(common_id, common_id) if current_gene_search_api else common_id
                    # Use primary_id in label to distinguish multiple ENS for same IGDB
                    expanded_labels.append(f"{gene_name} ({pid})")

            if expanded_rows:
                matrix = pd.DataFrame(expanded_rows)
                matrix.index = expanded_labels

            original_gene_ids = target_gene_ids  # Keep primary_ids for reference
            logger.info(f"  📊 Expanded matrix: {matrix.shape}")

            # Create visualization based on type
            transform_titles = {
                'raw': 'Raw TPM Values',
                'log2': 'Log2(TPM + 1)',
                'zscore_rows': 'Z-score (by gene)'
            }

            if viz_type == 'heatmap':
                # Create heatmap
                fig = go.Figure(data=go.Heatmap(
                    z=matrix.values,
                    x=[col.split('::')[1] for col in matrix.columns],  # Show sample names (format: Source::Sample)
                    y=list(matrix.index),
                    colorscale='Viridis',
                    showscale=True,
                    text=matrix.values.round(2) if show_values else None,
                    texttemplate="%{text}" if show_values else None,
                    textfont={"size": 8},
                    hoverongaps=False,
                    hovertemplate='<b>Gene:</b> %{y}<br><b>Sample:</b> %{x}<br><b>Expression:</b> %{z:.2f}<extra></extra>'
                ))

                fig.update_layout(
                    title=f"Expression Heatmap - {transform_titles.get(transform_method, transform_method)}",
                    xaxis_title="Samples",
                    yaxis_title="Genes",
                    height=max(400, len(matrix.index) * 30 + 150),  # Dynamic height based on gene count
                    margin=dict(l=200, r=50, t=100, b=50),
                    xaxis=dict(tickangle=45),
                    font=dict(size=10)
                )

            else:  # lineplot
                # Prepare data for line plot
                df_long = matrix.reset_index().melt(id_vars='index', var_name='Sample', value_name='Expression')
                df_long['Gene'] = df_long['index']
                df_long['Sample_Clean'] = df_long['Sample'].str.split('::').str[1]

                # Create line plot
                fig = px.line(df_long, x='Sample_Clean', y='Expression', color='Gene',
                             title=f"Expression Line Plot - {transform_titles.get(transform_method, transform_method)}",
                             markers=True)

                fig.update_layout(
                    xaxis_title="Samples",
                    yaxis_title=f"Expression ({transform_titles.get(transform_method, transform_method)})",
                    height=500,
                    margin=dict(l=50, r=50, t=100, b=100),
                    xaxis=dict(tickangle=45),
                    legend=dict(orientation="v", yanchor="top", y=1, xanchor="left", x=1.02),
                    font=dict(size=10)
                )

            # Create info panel with markdown
            info_text = f"""
### Expression Data Summary

**Dataset Overview:**
- **Genes:** {metadata['genes_found']}/{metadata['genes_requested']} found
- **Samples:** {metadata['samples_count']}
- **Sources:** {len(selected_sources)}

**Expression Statistics:**
- **Range:** {metadata['tpm_range'][0]:.2f} - {metadata['tpm_range'][1]:.2f} TPM
- **Coverage:** {metadata['non_zero_values']:,}/{metadata['total_values']:,} non-zero values ({metadata['non_zero_values']/metadata['total_values']*100:.1f}%)
"""

            info_panel = dbc.Card([
                dbc.CardBody([
                    dcc.Markdown(info_text, className="markdown-container")
                ])
            ], className="mb-3")

            # Handle missing genes
            missing_genes_alert = None
            if metadata['missing_genes']:
                missing_genes_alert = dbc.Alert([
                    html.Strong("Missing genes: "),
                    ", ".join(metadata['missing_genes'])
                ], color="warning", className="mb-3")

            components = [info_panel]

            if missing_genes_alert:
                components.append(missing_genes_alert)
            components.append(dcc.Graph(figure=fig, style={"width": "100%"}))

            # Prepare matrix data for export (convert to JSON-serializable format)
            export_data = {
                'matrix': matrix.to_dict('index'),
                'metadata': metadata,
                'transform_method': transform_method,
                'selected_sources': selected_sources,
                'selected_samples': selected_samples
            }

            return html.Div(components), export_data

        except Exception as e:
            logger.error(f"  ❌ Error in bulk RNA processing: {e}")
            return dbc.Alert(f"Error processing expression data: {str(e)}", color="danger", className="mt-4"), {}

    @app.callback(
        Output("download-expression-data", "data"),
        Input("export-csv-btn", "n_clicks"),
        State("expression-matrix-store", "data"),
        prevent_initial_call=True
    )
    def export_expression_data(csv_clicks, matrix_data):
        if not csv_clicks or not matrix_data:
            return no_update

        try:
            # Reconstruct DataFrame from stored data
            matrix_dict = matrix_data.get('matrix', {})
            if not matrix_dict:
                return no_update

            df = pd.DataFrame.from_dict(matrix_dict, orient='index')

            # Add metadata as header rows
            metadata = matrix_data.get('metadata', {})
            transform_method = matrix_data.get('transform_method', 'unknown')

            # Create export filename with timestamp
            timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")

            # Create CSV with metadata header
            output = []
            output.append(f"# Expression Data Export")
            output.append(f"# Export Date: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
            output.append(f"# Transform Method: {transform_method}")
            output.append(f"# Genes Found: {metadata.get('genes_found', 0)}")
            output.append(f"# Samples: {metadata.get('samples_count', 0)}")
            output.append(f"# Sources: {', '.join(metadata.get('sources_included', []))}")
            output.append("")  # Empty line before data

            # Add the data
            csv_data = df.to_csv(index_label="Gene")
            output.append(csv_data)

            return dict(
                content="\n".join(output),
                filename=f"expression_data_{timestamp}.csv",
                type="text/csv"
            )

        except Exception as e:
            logger.error(f"❌ Error in export: {e}")
            return no_update
