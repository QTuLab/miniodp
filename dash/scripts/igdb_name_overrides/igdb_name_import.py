#!/usr/bin/env python3
"""
Import gene-name mappings into the IGDB_NAME table (overwrite mode).
"""
import sqlite3
import sys
import pandas as pd
from pathlib import Path


SCRIPTS_DIR = Path(__file__).resolve().parents[1]
if str(SCRIPTS_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPTS_DIR))

from logging_utils import configure_script_logger, create_print_proxy  # type: ignore

_SCRIPT_LOGGER = configure_script_logger(__name__)
print = create_print_proxy(_SCRIPT_LOGGER)


def validate_data(df):
    """Validate TSV content before importing."""
    print(f"📊 Validating {len(df)} records...")
    
    # Required columns
    required_cols = ['IGDB', 'GeneName']
    missing_cols = [col for col in required_cols if col not in df.columns]
    if missing_cols:
        print(f"❌ Missing required columns: {missing_cols}")
        return False
    
    # Duplicate IGDB IDs (show up to 5)
    duplicates = df[df['IGDB'].duplicated()]
    if not duplicates.empty:
        print(f"⚠️  Found {len(duplicates)} duplicate IGDBs:")
        for igdb in duplicates['IGDB'].values[:5]:  # limit preview
            print(f"   - {igdb}")
        if len(duplicates) > 5:
            print(f"   ... and {len(duplicates) - 5} more")
    
    # Empty gene names
    empty_names = df[df['GeneName'].str.strip() == '']
    if not empty_names.empty:
        print(f"⚠️  Found {len(empty_names)} empty gene names")
    
    # Summary
    unique_genes = df['GeneName'].nunique()
    print(f"📈 Statistics:")
    print(f"   - Total records: {len(df)}")
    print(f"   - Unique IGDBs: {df['IGDB'].nunique()}")
    print(f"   - Unique gene names: {unique_genes}")
    print(f"   - Duplicates: {len(duplicates)}")
    print(f"   - Empty names: {len(empty_names)}")
    
    return len(duplicates) == 0 and len(empty_names) == 0

def import_igdb_names():
    """Import gene-name pairs into IGDB_NAME and overwrite old data."""
    
    # File locations
    dash_dir = Path(__file__).resolve().parents[2]
    db_path = dash_dir / "data" / "medaka" / "geneinfo.db"
    tsv_path = dash_dir / "data" / "medaka" / "misc" / "igdb_names.tsv"
    
    # Ensure database exists
    if not db_path.exists():
        print(f"❌ Database not found: {db_path}")
        print("💡 Please run create_igdb_info.py first")
        return False
    
    # Ensure TSV exists
    if not tsv_path.exists():
        print(f"❌ File not found: {tsv_path}")
        print("💡 Please create the igdb_names.tsv file first")
        return False
    
    # Load TSV file
    try:
        df = pd.read_csv(tsv_path, sep='\t', dtype=str, na_filter=False)
        print(f"📖 Loaded {len(df)} records from {tsv_path}")
        print(f"📋 Columns: {list(df.columns)}")
    except Exception as e:
        print(f"❌ Error reading TSV: {e}")
        return False
    
    # Validate content
    if not validate_data(df):
        print("❌ Data validation failed")
        return False
    
    # Add optional columns with safe defaults
    if 'Description' not in df.columns:
        df['Description'] = ''
    
    # Keep required columns only
    df = df[['IGDB', 'GeneName', 'Description']].copy()
    
    # Write to database (overwrite existing records)
    try:
        with sqlite3.connect(db_path) as conn:
            # Table exists?
            tables = pd.read_sql_query("SELECT name FROM sqlite_master WHERE type='table' AND name='IGDB_NAME'", conn)
            if tables.empty:
                print("❌ IGDB_NAME table not found. Please run igdb_name_create.py first")
                return False
            
            # Capture current count
            old_count = conn.execute("SELECT COUNT(*) FROM IGDB_NAME").fetchone()[0]
            print(f"📊 Current records in database: {old_count}")
            
            # Truncate table
            conn.execute("DELETE FROM IGDB_NAME")
            print("🗑️  Cleared existing data")
            
            # Insert new records
            df.to_sql('IGDB_NAME', conn, if_exists='append', index=False)
            
            # Report new totals
            new_count = conn.execute("SELECT COUNT(*) FROM IGDB_NAME").fetchone()[0]
            print(f"✅ Successfully imported {new_count} gene names")
            print(f"📈 Change: {old_count} → {new_count} ({new_count - old_count:+d})")
            
            return True
            
    except Exception as e:
        print(f"❌ Error importing to database: {e}")
        return False

if __name__ == "__main__":
    success = import_igdb_names()
    if success:
        print("\n🎉 Import completed successfully!")
    else:
        print("\n💥 Import failed!")
        exit(1)
