import dash_bootstrap_components as dbc
from dash import dcc, html

from .settings import build_prefixed_href


def create_sidebar(tf_families_list: str, default_search_field_options):
    return html.Div(
        [
            # Standalone toggle button section
            html.Div(
                [
                    dbc.Button(
                        "«",
                        id="sidebar-toggle-btn",
                        className="bg-transparent border-0 text-body-secondary",
                        style={
                            "position": "absolute",
                            "top": "10px",
                            "right": "10px",
                            "zIndex": "1001",
                            "fontSize": "14px",
                            "width": "30px",
                            "height": "30px",
                        },
                    )
                ],
                className="position-relative",
            ),
            # Collapsible header area
            html.Div(
                [
                    html.H5("Control Panel", className="mb-0 text-center"),
                    html.Div(id="current-species-display", className="text-center mt-2"),
                ],
                id="sidebar-header",
                className="p-3 border-bottom",
                style={"display": "block"},
            ),
            html.Div(
                [
                    html.Div(
                        [
                            html.Label("Search by Field", className="form-label fw-bold mb-2"),
                            dbc.Select(
                                id="search-field-select",
                                options=default_search_field_options,
                                value=default_search_field_options[0]["value"],
                                className="mb-2",
                            ),
                            html.Div(
                                [
                                    dcc.Input(
                                        id="search-terms-input",
                                        placeholder="Enter search term(s)...",
                                        className="form-control",
                                    ),
                                    dbc.Button(
                                        "View Options",
                                        id="view-options-btn",
                                        color="secondary",
                                        size="sm",
                                        className="mt-1",
                                        style={"display": "none"},
                                    ),
                                ],
                                className="mb-2",
                            ),
                            html.Div(id="search-results-container", className="mt-3"),
                        ],
                        id="gene-search-section",
                    ),
                    html.Div(
                        [
                            html.Div(
                                [
                                    html.Label("Target Genes", className="form-label fw-bold"),
                                    dbc.Button(
                                        "Clear All",
                                        id={"type": "clear-all-btn", "index": "clear"},
                                        color="danger",
                                        outline=True,
                                        size="sm",
                                        className="ms-3",
                                    ),
                                ],
                                className="d-flex align-items-center mb-2",
                            ),
                            html.Div(
                                id="gene-set-panel",
                                children=[
                                    dbc.Alert(
                                        "No genes in target list",
                                        color="info",
                                        className="mb-0 small",
                                    )
                                ],
                            ),
                        ],
                        className="mb-4 mt-3",
                    ),
                    html.Hr(),
                    html.Div(
                        id="contextual-controls-container",
                        style={
                            "maxHeight": "60vh",  # Cap the list height at 60% of viewport
                            "overflowY": "auto",  # Enable vertical scroll
                            "paddingRight": "5px",  # Leave room for the scrollbar
                        },
                    ),
                ],
                id="sidebar-content",
                className="p-3",
                style={
                    "display": "block",  # Visible by default
                    "overflowY": "auto",
                    "height": "calc(100vh - 120px)",  # Leave some breathing room
                    "maxHeight": "calc(100vh - 120px)",  # Match the calculated height
                    "paddingBottom": "20px",  # Reserve space at the bottom
                },
            ),
            dbc.Modal(
                [
                    dbc.ModalHeader(dbc.ModalTitle("Transcription Factor Families")),
                    dbc.ModalBody(
                        [
                            dcc.Markdown(
                                tf_families_list,
                                className="tf-families-container markdown-container",
                            )
                        ]
                    ),
                ],
                id="tf-modal",
                is_open=False,
                size="lg",
            ),
        ],
        id="sidebar",
        className="border-end",
        style={
            "width": "350px",
            "height": "100vh",
            "position": "fixed",
            "left": "0",
            "top": "0",
            "zIndex": "1000",
            "transition": "width 0.3s ease",
        },
    )


def create_main_content():
    return html.Div(
        [
            dcc.Tabs(
                id="main-tabs",
                value="tab-summary",
                className="custom-tabs",
                children=[
                    dcc.Tab(
                        label=name,
                        value=val,
                        className="custom-tab",
                        selected_className="custom-tab-selected",
                    )
                    for name, val in [
                        ("GeneInfo", "tab-summary"),
                        ("BulkRNA", "tab-bulk-rna"),
                        ("scRNA", "tab-scrna"),
                        ("scATAC", "tab-scatac"),
                        ("BulkMulti", "tab-bulkmulti"),
                    ]
                ],
            ),
            html.Div(id="tab-content", className="p-4"),
        ],
        id="main-content",
        style={"marginLeft": "350px", "minHeight": "100vh", "transition": "margin-left 0.3s ease"},
    )


def create_tab_content(tab_id):
    if tab_id == "tab-summary":
        return html.Div([dbc.Spinner(html.Div(id="gene-summary-content"))])
    if tab_id == "tab-bulk-rna":
        return html.Div([dbc.Spinner(html.Div(id="bulk-rna-content"))])
    if tab_id == "tab-scrna":
        return html.Div([dbc.Spinner(html.Div(id="scrna-content"))])
    if tab_id == "tab-scatac":
        return html.Div([dbc.Spinner(html.Div(id="scatac-content"))])
    if tab_id == "tab-bulkmulti":
        return html.Div(
            [
                dcc.Loading(
                    id="bulkmulti-loading",
                    type="circle",
                    delay_show=0,
                    children=html.Div(
                        id="peak2gene-content",
                        className="p2g-loading-container",
                        style={"minHeight": "420px"},
                    ),
                    fullscreen=False,
                )
            ],
            className="position-relative",
        )
    return html.Div([html.H3(f"Content for {tab_id}", className="text-primary mb-4")])


def create_species_selection_layout(available_species: list[str]):
    """Create species selection layout for root path."""
    return html.Div(
        [
            dbc.Container(
                [
                    html.H2("Multi-Omics Data Portal", className="text-center mb-4"),
                    html.P("Please select a species to explore:", className="text-center mb-4"),
                    dbc.Row(
                        [
                            dbc.Col(
                                [
                                    dbc.Card(
                                        [
                                            dbc.CardBody(
                                                [
                                                    html.H5(
                                                        species.replace("_", " ").title(),
                                                        className="card-title",
                                                    ),
                                                    html.P(
                                                        f"Explore {species} genomic data",
                                                        className="card-text",
                                                    ),
                                                    dbc.Button(
                                                        "View Data",
                                                        href=build_prefixed_href(species),
                                                        color="primary",
                                                        external_link=False,
                                                    ),
                                                ]
                                            )
                                        ],
                                        className="mb-3",
                                    )
                                ],
                                width=6,
                            )
                            for species in available_species
                        ],
                        justify="center",
                    ),
                    html.Hr(),
                    html.P(
                        f"Available species: {len(available_species)}",
                        className="text-center text-muted",
                    ),
                ],
                className="mt-5",
            )
        ]
    )


def create_species_error_layout(species_key: str, available_species: list[str]):
    """Create error layout for unavailable species."""
    return html.Div(
        [
            dbc.Container(
                [
                    dbc.Alert(
                        [
                            html.H4("Species Not Found", className="alert-heading"),
                            html.P(f"The species '{species_key}' is not available or data is not loaded."),
                            html.Hr(),
                            html.P(f"Available species: {', '.join(available_species)}", className="mb-0"),
                            html.Hr(),
                            dbc.Button(
                                "Back to Species Selection",
                                href=build_prefixed_href(),
                                color="primary",
                                className="mt-2",
                            ),
                        ],
                        color="warning",
                    )
                ],
                className="mt-5",
            )
        ]
    )


def create_layout():
    return html.Div(
        [
            # URL routing component
            dcc.Location(id="url", refresh=False),
            # Storage components
            dcc.Store(id="sidebar-state", data={"collapsed": False}),
            dcc.Store(id="target-gene-list-store", data=[]),
            dcc.Store(id="debounced-search-store"),
            dcc.Store(id="expression-matrix-store", data={}),
            dcc.Store(id="current-species-store", data=None),
            dcc.Store(id="scrna-page-store", data=0),
            dcc.Store(id="scatac-page-store", data=0),
            dcc.Store(id="bulkmulti-render-request"),
            # Page content (populated by routing callback)
            html.Div(
                id="page-content",
                children=[
                    # Default loading content
                    html.Div(
                        [
                            dbc.Spinner(
                                [
                                    html.H4("Multi-Omics Data Portal", className="text-center mt-5"),
                                    html.P(
                                        "Please navigate to a specific species URL",
                                        className="text-center text-muted",
                                    ),
                                ],
                                size="lg",
                            )
                        ],
                        className="d-flex justify-content-center align-items-center",
                        style={"height": "100vh"},
                    )
                ],
            ),
        ]
    )

