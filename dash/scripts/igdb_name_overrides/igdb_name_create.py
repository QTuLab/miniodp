#!/usr/bin/env python3
"""
Create the IGDB_NAME table (one-time utility).
"""
import sqlite3
import sys
from pathlib import Path


SCRIPTS_DIR = Path(__file__).resolve().parents[1]
if str(SCRIPTS_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPTS_DIR))

from logging_utils import configure_script_logger, create_print_proxy  # type: ignore

_SCRIPT_LOGGER = configure_script_logger(__name__)
print = create_print_proxy(_SCRIPT_LOGGER)

def create_igdb_name_table():
    """Create the IGDB_NAME lookup table if it is missing."""
    dash_dir = Path(__file__).resolve().parents[2]
    db_path = dash_dir / "data" / "medaka" / "geneinfo.db"
    
    if not db_path.exists():
        print(f"❌ Database not found: {db_path}")
        return False
    
    with sqlite3.connect(db_path) as conn:
        # Create table
        conn.execute("""
            CREATE TABLE IF NOT EXISTS IGDB_NAME (
                IGDB TEXT PRIMARY KEY,
                GeneName TEXT NOT NULL,
                Description TEXT
            )
        """)
        
        # Create supporting index
        conn.execute("CREATE INDEX IF NOT EXISTS idx_igdb_name_genename ON IGDB_NAME(GeneName)")
        
        conn.commit()
        print("✅ IGDB_NAME table created successfully")
        
        # Report current record count
        cursor = conn.execute("SELECT COUNT(*) FROM IGDB_NAME")
        count = cursor.fetchone()[0]
        print(f"📊 Current records in IGDB_NAME: {count}")
        
        return True

if __name__ == "__main__":
    create_igdb_name_table()
