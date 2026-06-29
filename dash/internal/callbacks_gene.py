from __future__ import annotations

import json
import logging

import dash
import dash_bootstrap_components as dbc
from dash import Input, Output, State, callback_context, html, no_update

from .apis import get_gene_search_api


def register(app, *, logger: logging.Logger) -> None:
    """Register gene search + target list callbacks."""

    @app.callback(
        Output("search-results-container", "children"),
        Input("debounced-search-store", "data"),
        State("search-field-select", "value"),
        State("current-species-store", "data"),
        prevent_initial_call=True,
    )
    def update_search_results(terms, field, species_key):
        if not terms or not terms.strip():
            return dbc.Alert("Enter a search term to see results.", color="info")

        current_gene_search_api = get_gene_search_api(species_key)
        if current_gene_search_api is None:
            return dbc.Alert("Search API is not available.", color="danger")

        results = current_gene_search_api.search_genes(terms, field=field, limit=100)
        logger.info(f"🔍 Search results for '{terms}' in field '{field}' on species '{species_key}':")
        for i, result in enumerate(results[:3]):  # Only show first 3 for brevity
            logger.info(
                f"  Result {i+1}: primary_id={result.get('primary_id')}, gene_name={result.get('gene_name')}"
            )

        if not results:
            return dbc.Alert("No results found.", color="info")

        all_result_ids = [res.get("primary_id", "") for res in results]

        result_items = []
        current_gene_search_api = get_gene_search_api(species_key)

        for res in results:
            primary_id = res.get("primary_id", "")

            if current_gene_search_api:
                display_config = current_gene_search_api.adapter.format_search_result_display(res)
                main_display = html.Strong(display_config["main_text"], className="d-block")

                secondary_info = []
                for item in display_config.get("secondary_items", []):
                    if item:
                        secondary_info.append(html.Small(item["text"], className=item["className"]))
            else:
                gene_name = res.get("gene_name", "N/A")
                primary_id_display = res.get("primary_id", "")
                main_display = html.Strong(gene_name, className="d-block")
                secondary_info = (
                    [html.Small(primary_id_display, className="text-muted")] if primary_id_display else []
                )

            result_item = dbc.ListGroupItem(
                [
                    dbc.Row(
                        [
                            dbc.Col([main_display] + secondary_info, width=9),
                            dbc.Col(
                                dbc.Button(
                                    "+",
                                    id={"type": "add-gene-btn", "index": primary_id},
                                    color="success",
                                    outline=True,
                                    size="sm",
                                ),
                                width=3,
                                className="d-flex align-items-center justify-content-end",
                            ),
                        ]
                    )
                ]
            )
            result_items.append(result_item)

        add_all_button = dbc.Button(
            "Add All",
            id={"type": "add-all-btn", "index": json.dumps(all_result_ids)},
            color="info",
            outline=True,
            size="sm",
            className="w-100 mb-2",
        )

        return [add_all_button, dbc.ListGroup(result_items, style={"maxHeight": "300px", "overflowY": "auto"})]

    @app.callback(
        Output("target-gene-list-store", "data"),
        Output("gene-set-panel", "children"),
        Input({"type": "add-gene-btn", "index": dash.ALL}, "n_clicks"),
        Input({"type": "remove-gene-btn", "index": dash.ALL}, "n_clicks"),
        Input({"type": "add-all-btn", "index": dash.ALL}, "n_clicks"),
        Input({"type": "clear-all-btn", "index": dash.ALL}, "n_clicks"),
        Input("current-species-store", "data"),
        State("target-gene-list-store", "data"),
        prevent_initial_call=True,
    )
    def manage_target_gene_list(
        add_clicks, remove_clicks, add_all_clicks, clear_all_clicks, species_key, current_target_list
    ):
        ctx = callback_context
        if not ctx.triggered:
            return no_update

        triggered_entry = ctx.triggered[0]
        triggered_prop_id = triggered_entry["prop_id"]
        if triggered_prop_id == "current-species-store.data":
            return [], dbc.Alert(
                "Species changed, target gene list cleared.", color="info", className="mb-0 small"
            )

        if not any(
            c
            for c in (add_clicks or []) + (remove_clicks or []) + (add_all_clicks or []) + (clear_all_clicks or [])
            if c is not None
        ):
            return no_update

        try:
            component_id = json.loads(triggered_prop_id.split(".")[0])
            component_type = component_id["type"]
        except (ValueError, KeyError):
            return no_update

        triggered_value = triggered_entry.get("value")
        if component_type != "clear-all-btn" and not triggered_value:
            return no_update

        if component_type == "clear-all-btn":
            return [], dbc.Alert("No genes in target list", color="info", className="mb-0 small")

        target_list = list(current_target_list or [])
        gene_id = component_id["index"]

        if component_type == "add-gene-btn":
            if gene_id not in target_list:
                target_list.append(gene_id)
        elif component_type == "remove-gene-btn":
            if gene_id in target_list:
                target_list.remove(gene_id)  # remove a single occurrence
        elif component_type == "add-all-btn":
            new_genes = json.loads(gene_id)
            for g in new_genes:
                if g not in target_list:
                    target_list.append(g)

        if not target_list:
            return [], dbc.Alert("No genes in target list", color="info", className="mb-0 small")

        current_gene_search_api = get_gene_search_api(species_key)
        id_to_name_map = current_gene_search_api.get_gene_names_by_ids(target_list) if current_gene_search_api else {}
        target_list.sort(key=lambda gid: (id_to_name_map.get(gid, gid) or gid).lower())

        rows = []
        for i in range(0, len(target_list), 2):
            chunk = target_list[i : i + 2]
            cols = []
            for gene in chunk:
                badge = dbc.Badge(
                    [
                        f"{id_to_name_map.get(gene, gene)}",
                        dbc.Button(
                            "×",
                            id={"type": "remove-gene-btn", "index": gene},
                            color="light",
                            size="sm",
                            className="ms-2 border-0 bg-transparent text-white",
                        ),
                    ],
                    color="primary",
                    className="w-100 p-2 d-inline-flex align-items-center justify-content-between",
                    title=gene,
                )
                cols.append(dbc.Col(badge, width=6, className="mb-1"))
            rows.append(dbc.Row(cols, className="g-1 mb-1"))

        return target_list, html.Div(rows)

