from __future__ import annotations

import logging

from dash import Input, Output, State, no_update

from .apis import get_gene_search_api


def register(app, *, logger: logging.Logger) -> None:
    """Register sidebar/UI callbacks."""

    app.clientside_callback(
        """
        function(value) {
            if (!value) { return dash.no_update; }
            const a = new Promise(resolve => { setTimeout(() => { resolve(value); }, 350); });
            return a;
        }
        """,
        Output("debounced-search-store", "data"),
        Input("search-terms-input", "value"),
    )

    @app.callback(
        Output("sidebar", "style"),
        Output("main-content", "style"),
        Output("sidebar-content", "style"),
        Output("sidebar-header", "style"),
        Output("sidebar-toggle-btn", "children"),
        Output("sidebar-state", "data"),
        Input("sidebar-toggle-btn", "n_clicks"),
        State("sidebar-state", "data"),
    )
    def toggle_sidebar(n_clicks, sidebar_state):
        if n_clicks is None:
            collapsed = False
        else:
            collapsed = not (sidebar_state or {}).get("collapsed", False)

        if collapsed:
            sidebar_style = {"width": "60px"}
            main_style = {"marginLeft": "60px"}
            content_style = {"display": "none"}
            header_style = {"display": "none"}
            button_icon = "»"
        else:
            sidebar_style = {"width": "350px"}
            main_style = {"marginLeft": "350px"}
            content_style = {
                "display": "block",
                "overflowY": "auto",
                "height": "calc(100vh - 120px)",
                "maxHeight": "calc(100vh - 120px)",
                "paddingBottom": "20px",
            }
            header_style = {"display": "block"}
            button_icon = "«"

        sidebar_style.update(
            {
                "height": "100vh",
                "position": "fixed",
                "left": "0",
                "top": "0",
                "zIndex": "1000",
                "transition": "width 0.3s ease",
            }
        )
        main_style.update({"minHeight": "100vh", "transition": "margin-left 0.3s ease"})
        return sidebar_style, main_style, content_style, header_style, button_icon, {"collapsed": collapsed}

    @app.callback(
        Output("view-options-btn", "style"),
        Input("search-field-select", "value"),
        State("current-species-store", "data"),
    )
    def toggle_view_options_button(field, species_key):
        """Check if view options should be displayed for the given search field."""
        try:
            current_gene_search_api = get_gene_search_api(species_key)
            if current_gene_search_api:
                adapter = current_gene_search_api.adapter
                return {"display": "block"} if adapter.is_tf_search_field(field) else {"display": "none"}
        except Exception as exc:
            logger.error(f"Error checking TF search field: {exc}")
        return {"display": "none"}

    @app.callback(
        Output("tf-modal", "is_open"),
        Input("view-options-btn", "n_clicks"),
        State("tf-modal", "is_open"),
        prevent_initial_call=True,
    )
    def toggle_tf_modal(n, is_open):
        return (not is_open) if n else is_open

