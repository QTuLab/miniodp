#!/usr/bin/env python3
"""
mini Omics Data Portal - Dash Application
"""

import argparse
import os

import dash
import dash_bootstrap_components as dbc
from dash import Input, Output

from internal.callbacks_bulk_rna import register as register_bulk_rna_callbacks
from internal.callbacks_bulkmulti import register as register_bulkmulti_callbacks
from internal.callbacks_contextual_controls import register as register_contextual_controls_callbacks
from internal.callbacks_gene import register as register_gene_callbacks
from internal.callbacks_gene_summary import register as register_gene_summary_callbacks
from internal.callbacks_routing import register as register_routing_callbacks
from internal.callbacks_scatac import register as register_scatac_callbacks
from internal.callbacks_scrna import register as register_scrna_callbacks
from internal.callbacks_sidebar import register as register_sidebar_callbacks
from internal.layout import (
    create_layout,
    create_tab_content,
)
from internal.logging import configure_logging
from internal.species_data import get_available_species
from internal.settings import DASH_DATA_PATH, DASH_PREFIX, load_tf_families_text


logger = configure_logging()

logger.info("🔧 Configuration:")
logger.info(f"  - DASH_PREFIX: {DASH_PREFIX}")
logger.info(f"  - DASH_DATA_PATH: {DASH_DATA_PATH}")

# --- App Initialization and Data Loading ---
app = dash.Dash(
    __name__,
    external_stylesheets=[dbc.themes.BOOTSTRAP],
    suppress_callback_exceptions=True,
    url_base_pathname=DASH_PREFIX
)
app.title = "mini Omics Data Portal"

server = app.server

TF_FAMILIES_LIST = load_tf_families_text()

AVAILABLE_SPECIES = get_available_species()
logger.info(f"📋 Available species: {AVAILABLE_SPECIES}")

DEFAULT_SEARCH_FIELD_OPTIONS = [{'label': 'Gene Name', 'value': 'gene_name'}]

app.layout = create_layout()

register_routing_callbacks(
    app,
    logger=logger,
    available_species=AVAILABLE_SPECIES,
    tf_families_list=TF_FAMILIES_LIST,
    default_search_field_options=DEFAULT_SEARCH_FIELD_OPTIONS,
)
register_sidebar_callbacks(app, logger=logger)
register_gene_callbacks(app, logger=logger)
register_gene_summary_callbacks(app, logger=logger)
register_contextual_controls_callbacks(app, logger=logger)
register_bulk_rna_callbacks(app, logger=logger)
register_bulkmulti_callbacks(app, logger=logger)
register_scrna_callbacks(app, logger=logger)
register_scatac_callbacks(app, logger=logger)

@app.callback(Output("tab-content", "children"), Input("main-tabs", "value"))
def update_tab_content(active_tab):
    logger.info(f"🔄 update_tab_content called with: {active_tab}")
    return create_tab_content(active_tab)

# --- Main Execution ---
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Start miniODP Gene Browser")
    parser.add_argument("--host", default=os.environ.get('DASH_HOST', '127.0.0.1'), help="Host address (use 0.0.0.0 to allow external access)")
    parser.add_argument("--port", type=int, default=int(os.environ.get('DASH_PORT', '8050')), help="Port number")
    parser.add_argument("--debug", action="store_true", help="Enable debug mode")
    args = parser.parse_args()
    
    # Override with environment variables if set
    host = os.environ.get('DASH_HOST', args.host)
    port = int(os.environ.get('DASH_PORT', args.port))
    debug = os.environ.get('DASH_DEBUG', 'false').lower() == 'true' or args.debug
    
    logger.info(f"🚀 Starting miniODP Gene Browser...")
    logger.info(f"🌍 Host: {host}")
    logger.info(f"🔌 Port: {port}")
    logger.info(f"🔧 Debug mode: {'Enabled' if debug else 'Disabled'}")
    app.run(host=host, port=port, debug=debug)
