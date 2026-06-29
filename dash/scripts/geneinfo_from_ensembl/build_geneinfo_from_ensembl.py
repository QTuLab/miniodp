#!/usr/bin/env python3
"""
Automated geneinfo.db builder based on Ensembl.
Usage: python build_geneinfo_from_ensembl.py danio_rerio

Author: miniodp project
Date: 2025-07-15
"""

import argparse
import sqlite3
import requests
import pandas as pd
import time
import sys
from io import StringIO
from tqdm import tqdm
from pathlib import Path
import logging

# Configure logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

class EnsemblGeneinfoBuilder:
    """Main helper that builds geneinfo.db from Ensembl."""
    
    def __init__(self, species, output_dir="."):
        self.species = species
        self.output_dir = Path(output_dir)
        
        # Output paths
        self.geneinfo_db = self.output_dir / "geneinfo.db"
        self.geneinfo_toml = self.output_dir / "geneinfo.toml"
        
        # API endpoints
        self.biomart_url = "https://www.ensembl.org/biomart/martservice"
        self.rest_url = "https://rest.ensembl.org"
        
        # Delay between requests to avoid hammering the APIs
        self.request_delay = 0.1
        
        # Ensure output directory exists
        self.output_dir.mkdir(parents=True, exist_ok=True)
    
    def build_biomart_query(self, dataset, attributes, filters=None):
        """Build BioMart query XML payload."""
        attrs = ''.join([f'<Attribute name="{attr}" />' for attr in attributes])
        
        filter_xml = ''
        if filters:
            filter_xml = ''.join([
                f'<Filter name="{name}" value="{value}" />' 
                for name, value in filters.items()
            ])
        
        return f"""<?xml version="1.0" encoding="UTF-8"?>
        <Query virtualSchemaName="default" formatter="TSV" header="1" uniqueRows="1">
            <Dataset name="{dataset}" interface="default">
                {attrs}
                {filter_xml}
            </Dataset>
        </Query>"""
    
    def fetch_biomart_data(self, dataset, attributes, filters=None, max_retries=3):
        """Fetch BioMart data with retry support."""
        query = self.build_biomart_query(dataset, attributes, filters)
        
        for attempt in range(max_retries):
            try:
                logger.info(f"Fetching data from {dataset} (attempt {attempt + 1}/{max_retries})")
                response = requests.post(
                    self.biomart_url, 
                    data={"query": query},
                    timeout=600  # 10-minute timeout
                )
                
                if response.status_code == 200:
                    # Check for inline error message
                    content = response.text.strip()
                    if content.startswith("Query ERROR") or len(content) == 0:
                        raise Exception(f"BioMart returned error or empty result: {content[:200]}")
                    
                    df = pd.read_csv(StringIO(content), sep='\t')
                    logger.info(f"Successfully fetched {len(df)} rows")
                    return df
                else:
                    raise Exception(f"HTTP {response.status_code}: {response.text[:200]}")
                    
            except Exception as e:
                logger.warning(f"Attempt {attempt + 1} failed: {e}")
                if attempt < max_retries - 1:
                    time.sleep(2 ** attempt)  # Exponential backoff
                else:
                    raise Exception(f"All {max_retries} attempts failed. Last error: {e}")
    
    def verify_species(self):
        """Verify that the requested species exists in Ensembl."""
        logger.info(f"Verifying species: {self.species}")
        
        try:
            # Query the species list once via REST
            url = f"{self.rest_url}/info/species"
            response = requests.get(url, params={"content-type": "application/json"}, timeout=30)
            
            if response.status_code == 200:
                species_data = response.json()
                species_list = [sp.get('name', '').lower() for sp in species_data.get('species', [])]
                
                if self.species.lower() in species_list:
                    logger.info(f"✅ Species {self.species} found in Ensembl")
                    return True
                else:
                    logger.error(f"❌ Species {self.species} not found in Ensembl")
                    logger.info(f"Available species: {', '.join(sorted(species_list)[:10])}...")
                    return False
            else:
                logger.warning("Could not verify species, proceeding anyway...")
                return True
                
        except Exception as e:
            logger.warning(f"Species verification failed: {e}, proceeding anyway...")
            return True
    
    def build_gene_info_table(self):
        """Build ENS_INFO table with gene metadata."""
        logger.info("📊 Fetching gene basic information...")
        
        dataset = f"{self.species}_gene_ensembl"
        attributes = [
            "ensembl_gene_id",
            "external_gene_name", 
            "gene_biotype",
            "description"
        ]
        
        df = self.fetch_biomart_data(dataset, attributes)
        
        # Basic cleanup
        df = df.dropna(subset=['ensembl_gene_id'])
        df['external_gene_name'] = df['external_gene_name'].fillna('')
        df['gene_biotype'] = df['gene_biotype'].fillna('unknown')
        df['description'] = df['description'].fillna('')
        
        # Normalize column names
        df = df.rename(columns={
            'ensembl_gene_id': 'ENS',
            'external_gene_name': 'GeneName',
            'gene_biotype': 'GeneType',
            'description': 'Description'
        })
        
        logger.info(f"✅ Gene info: {len(df)} genes")
        return df
    
    def build_locus_table(self):
        """Build ENS_Locus table with genomic coordinates."""
        logger.info("🧬 Fetching genomic coordinates...")
        
        dataset = f"{self.species}_gene_ensembl"
        attributes = [
            "ensembl_gene_id",
            "chromosome_name",
            "start_position",
            "end_position", 
            "strand"
        ]
        
        df = self.fetch_biomart_data(dataset, attributes)
        
        # Drop rows missing IDs
        df = df.dropna(subset=['ensembl_gene_id'])
        
        # Normalize chromosome names
        df['chromosome_name'] = df['chromosome_name'].astype(str)
        
        # Normalize column names
        df = df.rename(columns={
            'ensembl_gene_id': 'ENS',
            'chromosome_name': 'chromosome',
            'start_position': 'start',
            'end_position': 'end',
            'strand': 'strand'
        })
        
        logger.info(f"✅ Locus info: {len(df)} entries")
        return df
    
    def build_go_table(self):
        """Build ENS_GO table with GO annotations."""
        logger.info("🔬 Fetching GO annotations...")
        
        dataset = f"{self.species}_gene_ensembl"
        attributes = [
            "ensembl_gene_id",
            "go_id",
            "name_1006",  # GO term name
            "definition_1006",  # GO definition
            "go_linkage_type",  # Evidence code
            "namespace_1003"  # GO domain (BP/MF/CC)
        ]
        
        df = self.fetch_biomart_data(dataset, attributes)
        
        # Keep only valid GO accessions
        df = df.dropna(subset=['go_id'])
        df = df[df['go_id'].str.startswith('GO:', na=False)]
        
        # Fill missing text fields
        df['name_1006'] = df['name_1006'].fillna('')
        df['definition_1006'] = df['definition_1006'].fillna('')
        df['go_linkage_type'] = df['go_linkage_type'].fillna('IEA')
        df['namespace_1003'] = df['namespace_1003'].fillna('')
        
        # Normalize column names
        df = df.rename(columns={
            'ensembl_gene_id': 'ENS',
            'go_id': 'GO_Accession',
            'name_1006': 'GO_Name',
            'definition_1006': 'GO_Definition',
            'go_linkage_type': 'GO_Evidence',
            'namespace_1003': 'GO_Domain'
        })
        
        logger.info(f"✅ GO annotations: {len(df)} entries")
        return df
    
    def build_interpro_table(self):
        """Build ENS_INTERPRO table with InterPro domains."""
        logger.info("🧪 Fetching InterPro annotations...")
        
        dataset = f"{self.species}_gene_ensembl"
        attributes = [
            "ensembl_gene_id",
            "interpro",
            "interpro_description"
        ]
        
        try:
            df = self.fetch_biomart_data(dataset, attributes)
            
            # Drop invalid IDs
            df = df.dropna(subset=['interpro'])
            df = df[df['interpro'].str.startswith('IPR', na=False)]
            
            df['interpro_description'] = df['interpro_description'].fillna('')
            
            # Normalize column names
            df = df.rename(columns={
                'ensembl_gene_id': 'ENS',
                'interpro': 'IPRID',
                'interpro_description': 'IPRShortDesc'
            })
            
            logger.info(f"✅ InterPro annotations: {len(df)} entries")
            return df
            
        except Exception as e:
            logger.warning(f"InterPro data fetch failed: {e}")
            # Return empty frame with correct schema
            return pd.DataFrame(columns=['ENS', 'IPRID', 'IPRShortDesc'])
    
    def build_ortholog_table(self, target_species):
        """Build ortholog tables for selected reference species."""
        logger.info(f"🔗 Fetching {target_species} orthologs...")
        
        # Map helper names to Ensembl dataset codes
        species_map = {
            'human': 'hsapiens',
            'mouse': 'mmusculus', 
            'fly': 'dmelanogaster'
        }
        
        target_code = species_map.get(target_species, target_species)
        dataset = f"{self.species}_gene_ensembl"
        
        attributes = [
            "ensembl_gene_id",
            f"{target_code}_homolog_ensembl_gene",
            f"{target_code}_homolog_associated_gene_name",
            f"{target_code}_homolog_orthology_type"
        ]
        
        try:
            df = self.fetch_biomart_data(dataset, attributes)
            
            # Drop missing ortholog IDs
            ortho_col = f"{target_code}_homolog_ensembl_gene"
            df = df.dropna(subset=[ortho_col])
            df = df[df[ortho_col] != '']
            
            # Fill optional metadata
            name_col = f"{target_code}_homolog_associated_gene_name"
            type_col = f"{target_code}_homolog_orthology_type"
            
            df[name_col] = df[name_col].fillna('')
            df[type_col] = df[type_col].fillna('ortholog_one2one')
            
            # Normalize column names
            df = df.rename(columns={
                'ensembl_gene_id': 'ENS',
                ortho_col: 'OrthoID',
                name_col: 'OrthoName',
                type_col: 'OrthoType'
            })
            
            logger.info(f"✅ {target_species} orthologs: {len(df)} entries")
            return df
            
        except Exception as e:
            logger.warning(f"{target_species} ortholog data fetch failed: {e}")
            return pd.DataFrame(columns=['ENS', 'OrthoID', 'OrthoName', 'OrthoType'])
    
    def build_sequence_table(self):
        """Build ENS_Seq table with representative sequences."""
        logger.info("🧬 Fetching sequence information...")
        
        dataset = f"{self.species}_gene_ensembl"
        attributes = [
            "ensembl_gene_id",
            "ensembl_transcript_id",
            "gene_exon_intron"  # Genomic sequence data
        ]
        
        try:
            df = self.fetch_biomart_data(dataset, attributes)
            
            # Drop missing identifiers
            df = df.dropna(subset=['ensembl_gene_id', 'ensembl_transcript_id'])
            
            # Keep a single primary transcript per gene
            df = df.drop_duplicates(subset=['ensembl_gene_id'], keep='first')
            
            df['gene_exon_intron'] = df['gene_exon_intron'].fillna('')
            
            # Normalize column names
            df = df.rename(columns={
                'ensembl_gene_id': 'ENS',
                'ensembl_transcript_id': 'ENStrpt',
                'gene_exon_intron': 'nseq'
            })
            
            logger.info(f"✅ Sequence info: {len(df)} entries")
            return df
            
        except Exception as e:
            logger.warning(f"Sequence data fetch failed: {e}")
            return pd.DataFrame(columns=['ENS', 'ENStrpt', 'nseq'])
    
    def create_database_tables(self, conn):
        """Create database schema."""
        cursor = conn.cursor()
        
        # Drop pre-existing tables
        tables_to_drop = [
            'ENS_INFO', 'ENS_Locus', 'ENS_GO', 'ENS_INTERPRO', 
            'ENS_Seq', 'ENS_OrthoHuman', 'ENS_OrthoMouse', 'ENS_OrthoFly'
        ]
        
        for table in tables_to_drop:
            cursor.execute(f"DROP TABLE IF EXISTS {table}")
        
        # Recreate schema compatible with historical format
        table_schemas = {
            'ENS_INFO': '''
                CREATE TABLE ENS_INFO (
                    ENS TEXT,
                    GeneName TEXT,
                    GeneType TEXT,
                    Description TEXT
                )
            ''',
            'ENS_Locus': '''
                CREATE TABLE ENS_Locus (
                    ENS TEXT,
                    chromosome TEXT,
                    start INTEGER,
                    end INTEGER,
                    strand INTEGER
                )
            ''',
            'ENS_GO': '''
                CREATE TABLE ENS_GO (
                    ENS TEXT,
                    GO_Accession TEXT,
                    GO_Name TEXT,
                    GO_Definition TEXT,
                    GO_Evidence TEXT,
                    GO_Domain TEXT
                )
            ''',
            'ENS_INTERPRO': '''
                CREATE TABLE ENS_INTERPRO (
                    ENS TEXT,
                    IPRID TEXT,
                    IPRShortDesc TEXT
                )
            ''',
            'ENS_Seq': '''
                CREATE TABLE ENS_Seq (
                    ENS TEXT,
                    ENStrpt TEXT,
                    nseq TEXT
                )
            ''',
            'ENS_OrthoHuman': '''
                CREATE TABLE ENS_OrthoHuman (
                    ENS TEXT,
                    OrthoID TEXT,
                    OrthoName TEXT,
                    OrthoType TEXT
                )
            ''',
            'ENS_OrthoMouse': '''
                CREATE TABLE ENS_OrthoMouse (
                    ENS TEXT,
                    OrthoID TEXT,
                    OrthoName TEXT,
                    OrthoType TEXT
                )
            ''',
            'ENS_OrthoFly': '''
                CREATE TABLE ENS_OrthoFly (
                    ENS TEXT,
                    OrthoID TEXT,
                    OrthoName TEXT,
                    OrthoType TEXT
                )
            '''
        }
        
        for table_name, schema in table_schemas.items():
            cursor.execute(schema)
            logger.info(f"Created table: {table_name}")
        
        conn.commit()
    
    def create_indexes(self, conn):
        """Create database indexes to improve query performance."""
        logger.info("📈 Creating database indexes...")
        
        cursor = conn.cursor()
        
        indexes = [
            "CREATE INDEX IF NOT EXISTS idx_ens_info_ens ON ENS_INFO(ENS)",
            "CREATE INDEX IF NOT EXISTS idx_ens_info_genename ON ENS_INFO(GeneName)",
            "CREATE INDEX IF NOT EXISTS idx_ens_locus_ens ON ENS_Locus(ENS)",
            "CREATE INDEX IF NOT EXISTS idx_ens_go_ens ON ENS_GO(ENS)",
            "CREATE INDEX IF NOT EXISTS idx_ens_interpro_ens ON ENS_INTERPRO(ENS)",
            "CREATE INDEX IF NOT EXISTS idx_ens_seq_ens ON ENS_Seq(ENS)",
            "CREATE INDEX IF NOT EXISTS idx_ens_orthohuman_ens ON ENS_OrthoHuman(ENS)",
            "CREATE INDEX IF NOT EXISTS idx_ens_orthomouse_ens ON ENS_OrthoMouse(ENS)",
            "CREATE INDEX IF NOT EXISTS idx_ens_orthofly_ens ON ENS_OrthoFly(ENS)"
        ]
        
        created_count = 0
        for idx_sql in indexes:
            try:
                cursor.execute(idx_sql)
                created_count += 1
            except Exception as e:
                logger.warning(f"Index creation warning: {e}")
        
        conn.commit()
        logger.info(f"✅ Created {created_count} indexes")
    
    def insert_dataframe_to_table(self, conn, df, table_name):
        """Insert a DataFrame into the specified SQL table."""
        if len(df) > 0:
            df.to_sql(table_name, conn, if_exists='append', index=False)
            logger.info(f"✅ Inserted {len(df)} rows into {table_name}")
        else:
            logger.info(f"⚠️  No data to insert into {table_name}")
    
    def build_database(self):
        """Build the fully indexed geneinfo.db SQLite database."""
        logger.info(f"🚀 Building geneinfo.db for {self.species}")
        
        # Verify the requested species exists in the source metadata
        if not self.verify_species():
            return False
        
        # Drop any existing database file before rebuilding
        if self.geneinfo_db.exists():
            self.geneinfo_db.unlink()
            logger.info("Removed existing database file")
        
        try:
            with sqlite3.connect(self.geneinfo_db) as conn:
                # Create the schema first
                self.create_database_tables(conn)
                
                # 1. Core gene metadata
                gene_info = self.build_gene_info_table()
                self.insert_dataframe_to_table(conn, gene_info, 'ENS_INFO')
                
                if len(gene_info) == 0:
                    raise Exception("No gene information found - check species name")
                
                # 2. Genomic coordinates
                locus_info = self.build_locus_table()
                self.insert_dataframe_to_table(conn, locus_info, 'ENS_Locus')
                
                # 3. GO annotations
                go_info = self.build_go_table()
                self.insert_dataframe_to_table(conn, go_info, 'ENS_GO')
                
                # 4. InterPro domains
                interpro_info = self.build_interpro_table()
                self.insert_dataframe_to_table(conn, interpro_info, 'ENS_INTERPRO')
                
                # 5. Transcript sequences
                seq_info = self.build_sequence_table()
                self.insert_dataframe_to_table(conn, seq_info, 'ENS_Seq')
                
                # 6. Ortholog lookups
                ortholog_species = ['human', 'mouse', 'fly']
                for species in ortholog_species:
                    try:
                        ortho_info = self.build_ortholog_table(species)
                        table_name = f'ENS_Ortho{species.capitalize()}'
                        self.insert_dataframe_to_table(conn, ortho_info, table_name)
                    except Exception as e:
                        logger.warning(f"Failed to fetch {species} orthologs: {e}")
                
                # 7. Indexes
                self.create_indexes(conn)
                
                # 8. Refresh statistics
                cursor = conn.cursor()
                cursor.execute("ANALYZE")
                conn.commit()
                
                # 9. Summaries
                cursor.execute("SELECT COUNT(*) FROM ENS_INFO")
                gene_count = cursor.fetchone()[0]
                
                if gene_count > 0:
                    logger.info(f"✅ Database built successfully: {gene_count:,} genes")
                    
                    db_size = self.geneinfo_db.stat().st_size / (1024 * 1024)
                    logger.info(f"📁 Database size: {db_size:.1f} MB")
                    
                    return True
                
                raise Exception("Database built but contains no genes")
        except Exception as e:
            logger.error(f"❌ Error building database: {e}")
            return False
    
    def generate_toml_config(self):
        """Generate geneinfo.toml configuration."""
        logger.info("📝 Generating geneinfo.toml configuration...")
        
        # Base TOML template
        toml_content = [
            "# Geneinfo database search configuration",
            f"# Auto-generated for species: {self.species}",
            f"# Generated on: {time.strftime('%Y-%m-%d %H:%M:%S')}",
            "",
            "[[search_fields]]",
            'label = "Gene Name"',
            'value = "gene_name"',
            'db_table = "ENS_INFO"',
            'db_column = "GeneName"',
            'match_type = "like"',
            "",
            "[[search_fields]]",
            'label = "Ensembl ID"',
            'value = "ens_id"',
            'db_table = "ENS_INFO"',
            'db_column = "ENS"',
            'match_type = "exact"',
            "",
            "[[search_fields]]",
            'label = "Description"',
            'value = "description"',
            'db_table = "ENS_INFO"',
            'db_column = "Description"',
            'match_type = "like"',
            "",
            "[[search_fields]]",
            'label = "GO Term"',
            'value = "go_term"',
            'db_table = "ENS_GO"',
            'db_column = "GO_Name"',
            'match_type = "like"',
            "",
            "[[search_fields]]",
            'label = "GO Accession"',
            'value = "go_id"',
            'db_table = "ENS_GO"',
            'db_column = "GO_Accession"',
            'match_type = "exact"',
            "",
            "[[search_fields]]",
            'label = "Human Ortholog Name"',
            'value = "human_ortho_name"',
            'db_table = "ENS_OrthoHuman"',
            'db_column = "OrthoName"',
            'match_type = "like"',
            "",
            "[[search_fields]]",
            'label = "InterPro Domain"',
            'value = "interpro_desc"',
            'db_table = "ENS_INTERPRO"',
            'db_column = "IPRShortDesc"',
            'match_type = "like"',
        ]
        
        try:
            with open(self.geneinfo_toml, 'w', encoding='utf-8') as f:
                f.write('\n'.join(toml_content))
            
            logger.info(f"✅ Generated configuration: {self.geneinfo_toml}")
            return True
            
        except Exception as e:
            logger.error(f"❌ Error generating TOML config: {e}")
            return False
    
    def run_build(self):
        """Execute the full build workflow."""
        start_time = time.time()
        
        logger.info(f"🚀 Starting geneinfo.db build for {self.species}")
        logger.info(f"📂 Output directory: {self.output_dir}")
        
        try:
            # Build the SQLite database
            if not self.build_database():
                return False
            
            # Generate TOML config
            if not self.generate_toml_config():
                logger.warning("Database built successfully but TOML generation failed")
            
            # Report elapsed time
            elapsed = time.time() - start_time
            logger.info(f"🎉 Build completed successfully in {elapsed:.1f} seconds!")
            
            # List output files
            logger.info("📁 Output files:")
            if self.geneinfo_db.exists():
                db_size = self.geneinfo_db.stat().st_size / (1024 * 1024)
                logger.info(f"   - {self.geneinfo_db} ({db_size:.1f} MB)")
            
            if self.geneinfo_toml.exists():
                toml_size = self.geneinfo_toml.stat().st_size
                logger.info(f"   - {self.geneinfo_toml} ({toml_size} bytes)")
            
            return True
            
        except Exception as e:
            logger.error(f"❌ Build failed: {e}")
            return False

def main():
    """CLI entry point."""
    parser = argparse.ArgumentParser(
        description="Build geneinfo.db from Ensembl data",
        epilog="Example: python build_geneinfo_from_ensembl.py danio_rerio -o ./zebrafish_data"
    )
    
    parser.add_argument(
        "species", 
        help="Ensembl species name (e.g., danio_rerio, mus_musculus, homo_sapiens)"
    )
    parser.add_argument(
        "-o", "--output-dir", 
        default=".", 
        help="Output directory (default: current directory)"
    )
    parser.add_argument(
        "-v", "--verbose", 
        action="store_true", 
        help="Enable verbose logging"
    )
    
    args = parser.parse_args()
    
    # Adjust logging level
    if args.verbose:
        logging.getLogger().setLevel(logging.DEBUG)
    
    # Sanity-check required input
    if not args.species:
        logger.error("Species name is required")
        return 1
    
    # Instantiate builder and run
    builder = EnsemblGeneinfoBuilder(
        species=args.species,
        output_dir=args.output_dir,
    )
    
    success = builder.run_build()
    return 0 if success else 1

if __name__ == "__main__":
    sys.exit(main())
