"""Callbacks for gene summary page (detailed gene information display)."""

from __future__ import annotations

import json
import logging

import dash
import dash_bootstrap_components as dbc
from dash import Input, Output, State, callback_context, html, no_update

from .apis import get_gene_search_api
from .settings import build_jbrowse_url


EXTERNAL_LINK_PROPS = {"target": "_blank", "rel": "noopener noreferrer"}


def register(app, *, logger: logging.Logger) -> None:
    """Register gene summary callbacks."""

    @app.callback(
        Output("gene-summary-content", "children"),
        Input("target-gene-list-store", "data"),
        Input("summary-column-checklist", "value"),
        Input("main-tabs", "value"),
        State("current-species-store", "data")
    )
    def update_gene_summary_page(target_gene_ids, selected_columns, active_tab, species_key):
        logger.info(f"🔍 update_gene_summary_page called:")
        logger.info(f"  - active_tab: {active_tab}")
        logger.info(f"  - target_gene_ids: {target_gene_ids}")
        logger.info(f"  - selected_columns: {selected_columns}")

        if active_tab != "tab-summary":
            logger.error("  ❌ Not summary tab, returning no_update")
            return no_update

        if not target_gene_ids:
            logger.error("  ❌ No target genes")
            return dbc.Alert("Select genes from the sidebar to see their summary.", color="info", className="mt-4")

        current_gene_search_api = get_gene_search_api(species_key)
        if current_gene_search_api is None:
            logger.error("  ❌ API not available")
            return dbc.Alert("Gene search is not available for this species.", color="warning", className="mt-4")

        logger.info(f"  ✅ Calling get_gene_details for {len(target_gene_ids)} genes")
        details = current_gene_search_api.get_gene_details(target_gene_ids)
        logger.info(f"  📊 Retrieved details: {len(details)} records")
        logger.info(f"  📝 Sample data: {list(details.keys())[:3] if details else 'None'}")
        if details:
            first_gene = list(details.values())[0]
            logger.info(f"  🔍 First gene fields: {list(first_gene.keys())}")

        # JBrowse base URL from Hugo config
        jbrowse_url_base = None
        try:
            hugo_config = current_gene_search_api.adapter._load_hugo_species_config()
            jbrowse_url_base = hugo_config.get('jbrowse_url')
        except Exception as e:
            logger.warning(f"Failed to load JBrowse URL for {species_key}: {e}")

        if not details:
            logger.error("  ❌ No details retrieved")
            return dbc.Alert("Could not retrieve details for selected genes.", color="warning", className="mt-4")

        # 1. Create Simple Summary Table
        # Get column configuration from adapter instead of hardcoded values
        columns_config = current_gene_search_api.adapter.get_summary_columns()

        # Filter selected columns to only include valid ones for current species
        valid_columns = [col for col in selected_columns if col in columns_config]
        header = [html.Th(columns_config[col]) for col in valid_columns]

        rows = []
        for gene_id in target_gene_ids:
            gene_data = details.get(gene_id, {})

            # Truncate description for table display
            display_desc = gene_data.get('description') or 'N/A'
            if len(display_desc) > 100:
                display_desc = display_desc[:97] + '...'

            row_data = []
            for col in valid_columns:
                if col == 'description':
                    value = display_desc
                elif col == 'tf_families':
                    tf_families = gene_data.get('tf_families', [])
                    if tf_families:
                        # Create linked TF family badges for table
                        tf_links = []
                        for tf in tf_families:
                            tf_url = f"https://guolab.wchscu.cn/AnimalTFDB4/#/family_summary/TF/{tf}"
                            tf_link = html.A(
                                tf,
                                href=tf_url,
                                className="text-decoration-none me-1",
                                **EXTERNAL_LINK_PROPS,
                            )
                            tf_links.append(tf_link)
                            if tf != tf_families[-1]:  # Add comma if not last
                                tf_links.append(", ")
                        value = tf_links
                    else:
                        value = 'N/A'
                elif col == 'primary_id':
                    # Handle primary ID - check if it should have Ensembl link
                    primary_id = gene_data.get('primary_id') or 'N/A'
                    if primary_id != 'N/A':
                        # Check if primary_id is ENS format
                        if primary_id.startswith('ENS'):
                            external_links = gene_data.get('external_links', {})
                            ensembl_url = external_links.get('Ensembl')
                            if ensembl_url:
                                value = html.A(
                                    primary_id,
                                    href=ensembl_url,
                                    className="text-decoration-none",
                                    **EXTERNAL_LINK_PROPS,
                                )
                            else:
                                value = primary_id
                        else:
                            value = primary_id
                    else:
                        value = primary_id
                elif col == 'secondary_id':
                    # Handle secondary ID - check if it should have Ensembl link
                    secondary_id = gene_data.get('secondary_id') or 'N/A'
                    if secondary_id != 'N/A':
                        # Check if secondary_id is ENS format
                        if secondary_id.startswith('ENS'):
                            external_links = gene_data.get('external_links', {})
                            ensembl_url = external_links.get('Ensembl')
                            if ensembl_url:
                                value = html.A(
                                    secondary_id,
                                    href=ensembl_url,
                                    className="text-decoration-none",
                                    **EXTERNAL_LINK_PROPS,
                                )
                            else:
                                value = secondary_id
                        else:
                            value = secondary_id
                    else:
                        value = secondary_id
                elif col == 'locus':
                    locus_val = gene_data.get('locus') or 'N/A'
                    if locus_val != 'N/A' and jbrowse_url_base:
                        jb_url = build_jbrowse_url(jbrowse_url_base, locus_val)
                        value = html.A(
                            locus_val,
                            href=jb_url,
                            className="text-decoration-none",
                            **EXTERNAL_LINK_PROPS,
                        )
                    else:
                        value = locus_val
                else:
                    value = gene_data.get(col) or 'N/A'
                row_data.append(html.Td(value))
            rows.append(html.Tr(row_data))

        summary_table = dbc.Table([html.Thead(html.Tr(header)), html.Tbody(rows)], bordered=True, striped=True, hover=True, responsive=True, className="mb-5")

        # 2. Create Detail Cards
        cards = []
        for gene_id in target_gene_ids:
            gene_data = details.get(gene_id, {})

            # Accordion for details
            accordion_items = []
            if gene_data.get('go_terms'):
                # Filter out GO terms with None values
                valid_go_terms = [go for go in gene_data['go_terms'] if go.get('GO_Accession') and go.get('GO_Name')]
                if valid_go_terms:
                    go_table = dbc.Table(
                        [html.Thead(html.Tr([html.Th("Domain"), html.Th("ID"), html.Th("Name")]))] +
                        [html.Tbody([html.Tr([html.Td(go['GO_Domain'] or 'N/A'), html.Td(go['GO_Accession']), html.Td(go['GO_Name'])]) for go in valid_go_terms])],
                        size="sm", striped=True
                    )
                    accordion_items.append(dbc.AccordionItem(go_table, title="GO Annotation"))

            if gene_data.get('interpro'):
                # Remove duplicates from InterPro data
                seen_ipr = set()
                unique_interpro = []
                for ipr in gene_data['interpro']:
                    ipr_key = (ipr['IPRID'], ipr['IPRShortDesc'])
                    if ipr_key not in seen_ipr:
                        seen_ipr.add(ipr_key)
                        unique_interpro.append(ipr)

                if unique_interpro:
                    ipr_table = dbc.Table(
                        [html.Thead(html.Tr([html.Th("ID"), html.Th("Description")]))] +
                        [html.Tbody([html.Tr([html.Td(ipr['IPRID']), html.Td(ipr['IPRShortDesc'])]) for ipr in unique_interpro])],
                        size="sm", striped=True
                    )
                    accordion_items.append(dbc.AccordionItem(ipr_table, title="InterPro Domains"))

            # On-demand sequence loading
            accordion_items.append(dbc.AccordionItem(
                dbc.Spinner(html.Div(id={"type": "sequence-container", "index": gene_id}, children=[
                    dbc.Button("Load Sequence", id={"type": "load-seq-btn", "index": gene_id}, size="sm", color="light")
                ])),
                title="Sequence"
            ))

            # Build external database buttons dynamically from adapter-provided links
            external_links = gene_data.get('external_links', {})
            # JBrowse URL using Hugo configuration with locus padding
            if current_gene_search_api:
                hugo_config = current_gene_search_api.adapter._load_hugo_species_config()
                base_jbrowse_url = hugo_config.get('jbrowse_url', '')
                locus = gene_data.get('locus', '')
                jbrowse_url = build_jbrowse_url(base_jbrowse_url, locus) or base_jbrowse_url or '#'
            else:
                jbrowse_url = '#'

            header_buttons = []
            for link_name, link_url in external_links.items():
                if link_url and link_url != '#':
                    header_buttons.append(
                        dbc.Button(
                            link_name,
                            href=link_url,
                            color="primary",
                            outline=True,
                            size="sm",
                            external_link=True,
                            **EXTERNAL_LINK_PROPS,
                        )
                    )
            if jbrowse_url and jbrowse_url != '#':
                header_buttons.append(
                    dbc.Button(
                        "JBrowse",
                        href=jbrowse_url,
                        color="secondary",
                        outline=True,
                        size="sm",
                        external_link=True,
                        **EXTERNAL_LINK_PROPS,
                    )
                )

            card_header = html.Div([
                html.H5(gene_data.get('gene_name', gene_id), className="m-0"),
                dbc.ButtonGroup(header_buttons, className="ms-auto") if header_buttons else html.Div(className="ms-auto")
            ], className="d-flex align-items-center")

            # Build card body with basic info and TF families
            locus_value = gene_data.get('locus') or 'N/A'
            if locus_value != 'N/A' and jbrowse_url and jbrowse_url != '#':
                locus_component = html.A(
                    locus_value,
                    href=jbrowse_url,
                    className="text-decoration-none",
                    **EXTERNAL_LINK_PROPS,
                )
            else:
                locus_component = locus_value

            card_body_items = [
                dbc.Row([
                    dbc.Col(html.B("Locus:"), width="auto"), dbc.Col(locus_component),
                ], className="mb-1"),
                dbc.Row([
                    dbc.Col(html.B("Gene Type:"), width="auto"), dbc.Col(gene_data.get('gene_type') or 'N/A'),
                ], className="mb-1"),
            ]

            # Add TF families if present
            tf_families = gene_data.get('tf_families', [])
            if tf_families:
                # Create linked TF badges
                tf_badges = []
                for tf in tf_families:
                    tf_url = f"https://guolab.wchscu.cn/AnimalTFDB4/#/family_summary/TF/{tf}"
                    tf_badge = html.A(
                        dbc.Badge(tf, color="info", className="me-1"),
                        href=tf_url,
                        className="text-decoration-none",
                        **EXTERNAL_LINK_PROPS,
                    )
                    tf_badges.append(tf_badge)

                card_body_items.append(
                    dbc.Row([
                        dbc.Col(html.B("TF Families:"), width="auto"),
                        dbc.Col(tf_badges),
                    ], className="mb-2")
                )

            # Add ortholog information
            human_orthologs = gene_data.get('orthologs_human', [])
            valid_human = [orth for orth in human_orthologs if orth.get('OrthoID')]
            if valid_human:
                human_links = []
                for orth in valid_human:
                    link_text = orth.get('OrthoName') or orth['OrthoID']
                    human_links.append(
                        html.A(
                            link_text,
                            href=f"https://ensembl.org/Homo_sapiens/Gene/Summary?g={orth['OrthoID']}",
                            className="text-decoration-none me-2",
                            **EXTERNAL_LINK_PROPS,
                        )
                    )
                card_body_items.append(
                    dbc.Row([
                        dbc.Col(html.B("Human Orthologs:"), width="auto"),
                        dbc.Col(human_links),
                    ], className="mb-1")
                )

            mouse_orthologs = gene_data.get('orthologs_mouse', [])
            valid_mouse = [orth for orth in mouse_orthologs if orth.get('OrthoID')]
            if valid_mouse:
                mouse_links = []
                for orth in valid_mouse:
                    link_text = orth.get('OrthoName') or orth['OrthoID']
                    mouse_links.append(
                        html.A(
                            link_text,
                            href=f"https://ensembl.org/Mus_musculus/Gene/Summary?g={orth['OrthoID']}",
                            className="text-decoration-none me-2",
                            **EXTERNAL_LINK_PROPS,
                        )
                    )
                card_body_items.append(
                    dbc.Row([
                        dbc.Col(html.B("Mouse Orthologs:"), width="auto"),
                        dbc.Col(mouse_links),
                    ], className="mb-1")
                )

            card_body_items.extend([
                dbc.Card([
                    dbc.CardBody([
                        html.P([
                            html.B("Description: "),
                            gene_data.get('description', 'No description available.')
                        ])
                    ])
                ], className="mt-3"),
                dbc.Accordion(accordion_items, start_collapsed=True, flush=True, className="mt-3")
            ])

            card_body = card_body_items
            cards.append(dbc.Card([dbc.CardHeader(card_header), dbc.CardBody(card_body)], className="mb-3"))

        return html.Div([summary_table, html.Div(cards)])

    @app.callback(
        Output({"type": "sequence-container", "index": dash.ALL}, "children"),
        Input({"type": "load-seq-btn", "index": dash.ALL}, "n_clicks"),
        State("current-species-store", "data"),
        prevent_initial_call=True
    )
    def load_sequence(n_clicks, species_key):
        ctx = callback_context
        if not ctx.triggered or not any(n_clicks): return [no_update] * len(n_clicks)
        triggered_id_str = ctx.triggered[0]['prop_id'].split('.')[0]
        gene_id = json.loads(triggered_id_str)['index']
        current_gene_search_api = get_gene_search_api(species_key)
        sequences = current_gene_search_api.get_gene_sequence(gene_id) if current_gene_search_api else []

        outputs = [no_update] * len(ctx.outputs_list)
        triggered_index = -1
        for i, output in enumerate(ctx.outputs_list):
            if output['id']['index'] == gene_id:
                triggered_index = i
                break

        if triggered_index != -1:
            if not sequences:
                outputs[triggered_index] = html.Div("No sequences found.", className="text-muted")
            else:
                # Combine all sequences into one text block for easy copying
                all_sequences = "\n\n".join(seq['fasta'] for seq in sequences)

                # Add summary info if multiple sequences
                if len(sequences) > 1:
                    summary_text = f"Found {len(sequences)} transcript sequences:\n\n"
                    all_sequences = summary_text + all_sequences

                outputs[triggered_index] = html.Pre(
                    all_sequences,
                    style={
                        "maxHeight": "400px",
                        "overflowY": "auto",
                        "fontSize": "0.75em",
                        "fontFamily": "monospace",
                        "backgroundColor": "#f8f9fa",
                        "padding": "10px",
                        "border": "1px solid #dee2e6",
                        "borderRadius": "4px",
                        "whiteSpace": "pre-wrap"
                    }
                )

        return outputs
