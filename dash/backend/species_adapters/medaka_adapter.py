"""
Medaka species adapter for IGDB+ENS dual system.

This adapter handles medaka's unique data model where genes have both
IGDB and ENS identifiers, with expression data keyed by IGDB.
"""

import pandas as pd
import sqlite3
from pathlib import Path
from typing import List, Dict, Any
from .base_adapter import SpeciesAdapter


class MedakaAdapter(SpeciesAdapter):
    """Adapter for medaka with IGDB+ENS dual identifier system."""
    
    def get_search_fields(self) -> List[Dict[str, str]]:
        """Return search fields specific to medaka."""
        return [
            {'label': 'IGDB Gene ID', 'value': 'igdb_id'},
            {'label': 'Ensembl ID', 'value': 'ens_id'},
            {'label': 'Gene Name', 'value': 'gene_name'}, 
            {'label': 'GO ID', 'value': 'go_id'},
            {'label': 'IGDB TF Family', 'value': 'igdb_tf'},
            {'label': 'InterPro ID', 'value': 'interpro_id'},
            {'label': 'Human Ortholog ID', 'value': 'human_ortho_id'},
            {'label': 'Mouse Ortholog ID', 'value': 'mouse_ortho_id'},
            {'label': 'Zebrafish Ortholog ID', 'value': 'zebrafish_ortho_id'},
            {'label': 'Human Ortholog Name', 'value': 'human_ortho_name'},
            {'label': 'Mouse Ortholog Name', 'value': 'mouse_ortho_name'},
            {'label': 'Zebrafish Ortholog Name', 'value': 'zebrafish_ortho_name'},
        ]
    
    def search_genes(self, db_conn: sqlite3.Connection, query: str, search_type: str) -> pd.DataFrame:
        """
        Execute gene search for medaka with IGDB as primary identifier.
        
        This method handles both TOML configuration and adapter defaults,
        but always ensures IGDB is returned as the primary identifier.
        """
        if not query.strip():
            return pd.DataFrame()
        
        cleaned_query = self._clean_query(query)
        if not cleaned_query:
            return pd.DataFrame()
        
        # For gene_name searches, always use adapter method to search both ENS_INFO and IGDB_INFO
        if search_type == 'gene_name':
            return self._search_with_adapter_defaults(db_conn, cleaned_query, search_type)
        
        # Check if we have TOML configuration for this search type
        if search_type in self._search_field_configs:
            return self._search_with_toml_config(db_conn, cleaned_query, search_type)
        
        # Fallback to adapter defaults
        return self._search_with_adapter_defaults(db_conn, cleaned_query, search_type)
    
    def _search_with_toml_config(self, db_conn: sqlite3.Connection, query: str, search_type: str) -> pd.DataFrame:
        """Execute search using TOML configuration with IGDB priority."""
        config = self._search_field_configs[search_type]
        db_table = config['db_table']
        db_column = config['db_column']
        match_type = config.get('match_type', 'like')
        
        # For IGDB table searches
        if db_table.startswith("IGDB"):
            return self._search_igdb_table_with_toml(db_conn, query, db_table, db_column, match_type)
        
        # For ENS table searches (but return IGDB as primary)
        return self._search_ens_table_with_igdb_primary(db_conn, query, db_table, db_column, match_type)
    
    def _search_igdb_table_with_toml(self, db_conn: sqlite3.Connection, query: str, db_table: str, db_column: str, match_type: str) -> pd.DataFrame:
        """Search IGDB tables and return complete information."""
        base_select = """
        SELECT DISTINCT 
            l.IGDB, l.chromosome, l.start, l.end, l.strand,
            m.ENS, e.GeneName, e.GeneType, e.Description
        FROM IGDB_Locus l
        LEFT JOIN IGDB_ENS m ON l.IGDB = m.IGDB
        LEFT JOIN ENS_INFO e ON m.ENS = e.ENS
        """
        
        search_terms = [term.strip() for term in query.split() if term.strip()]
        all_results = pd.DataFrame()
        
        for term in search_terms:
            operator = "=" if match_type == "exact" else "LIKE"
            search_pattern = term if match_type == "exact" else f"%{term}%"
            
            if db_table == "IGDB_Locus":
                sql = f"{base_select} WHERE l.{db_column} {operator} ?"
            else:
                sql = f"{base_select} JOIN {db_table} t ON l.IGDB = t.IGDB WHERE t.{db_column} {operator} ?"
            
            sql += " ORDER BY l.IGDB"
            
            try:
                result_df = pd.read_sql_query(sql, db_conn, params=[search_pattern])
                if not result_df.empty:
                    all_results = pd.concat([all_results, result_df], ignore_index=True)
            except Exception as e:
                self.logger.warning("IGDB TOML search error for term '%s': %s", term, e)
        
        if not all_results.empty:
            return all_results.drop_duplicates().reset_index(drop=True)
        return pd.DataFrame()

    def _search_ens_table_with_igdb_primary(self, db_conn: sqlite3.Connection, query: str, db_table: str, db_column: str, match_type: str) -> pd.DataFrame:
        """Search ENS tables but return IGDB as primary for medaka."""
        base_select = """
        SELECT DISTINCT 
            m.IGDB, l.chromosome, l.start, l.end, l.strand,
            e.ENS, e.GeneName, e.GeneType, e.Description
        FROM ENS_INFO e
        JOIN IGDB_ENS m ON e.ENS = m.ENS
        JOIN IGDB_Locus l ON m.IGDB = l.IGDB
        """
        
        search_terms = [term.strip() for term in query.split() if term.strip()]
        all_results = pd.DataFrame()
        
        for term in search_terms:
            operator = "=" if match_type == "exact" else "LIKE"
            search_pattern = term if match_type == "exact" else f"%{term}%"
            
            if db_table == "ENS_INFO":
                sql = f"{base_select} WHERE e.{db_column} {operator} ?"
            else:
                sql = f"{base_select} JOIN {db_table} t ON e.ENS = t.ENS WHERE t.{db_column} {operator} ?"
            
            sql += " ORDER BY m.IGDB"
            
            try:
                result_df = pd.read_sql_query(sql, db_conn, params=[search_pattern])
                if not result_df.empty:
                    all_results = pd.concat([all_results, result_df], ignore_index=True)
            except Exception as e:
                self.logger.warning("ENS TOML search error for term '%s': %s", term, e)
        
        if not all_results.empty:
            return all_results.drop_duplicates().reset_index(drop=True)
        return pd.DataFrame()

    def batch_search_genes(self, db_conn: sqlite3.Connection, gene_ids: List[str], search_type: str) -> pd.DataFrame:
        """
        Execute gene search for multiple genes efficiently using IN clause for medaka species.
        
        This method is optimized for 'igdb_id' or 'ens_id' search types.
        For other types, it falls back to individual searches or raises an error.
        """
        if not gene_ids:
            return pd.DataFrame()

        placeholders = ','.join(['?' for _ in gene_ids])
        
        # Base select statement to retrieve all necessary information (IGDB as primary)
        base_select_igdb_primary = """
        SELECT DISTINCT 
            l.IGDB, l.chromosome, l.start, l.end, l.strand,
            m.ENS, e.GeneName, e.GeneType, e.Description
        FROM IGDB_Locus l
        LEFT JOIN IGDB_ENS m ON l.IGDB = m.IGDB
        LEFT JOIN ENS_INFO e ON m.ENS = e.ENS
        """
        
        if search_type == 'igdb_id' or search_type == 'primary_id':
            sql = f"{base_select_igdb_primary} WHERE l.IGDB IN ({placeholders}) ORDER BY l.IGDB"
            return pd.read_sql_query(sql, db_conn, params=gene_ids)
        elif search_type == 'ens_id' or search_type == 'secondary_id':
            # Need to search by ENS but return IGDB-primary results
            sql = f"""
            SELECT DISTINCT 
                m.IGDB, l.chromosome, l.start, l.end, l.strand,
                e.ENS, e.GeneName, e.GeneType, e.Description
            FROM ENS_INFO e
            JOIN IGDB_ENS m ON e.ENS = m.ENS
            JOIN IGDB_Locus l ON m.IGDB = l.IGDB
            WHERE e.ENS IN ({placeholders})
            ORDER BY m.IGDB
            """
            return pd.read_sql_query(sql, db_conn, params=gene_ids)
        elif search_type == 'gene_name':
            self.logger.warning(f"Batch search for 'gene_name' is not optimized; performing individual searches. Consider using ID-based search for performance.")
            all_results = pd.DataFrame()
            for gene_id in gene_ids: # gene_ids here are gene names to search
                result = self._search_by_gene_name(db_conn, gene_id)
                if not result.empty:
                    all_results = pd.concat([all_results, result], ignore_index=True)
            return all_results.drop_duplicates().reset_index(drop=True)
        else:
            self.logger.warning(f"Batch search type '{search_type}' not directly optimized for Medaka; performing individual searches.")
            all_results = pd.DataFrame()
            for gene_id in gene_ids:
                result = self.search_genes(db_conn, gene_id, search_type) # Calls the single-query search_genes
                if not result.empty:
                    all_results = pd.concat([all_results, result], ignore_index=True)
            return all_results.drop_duplicates().reset_index(drop=True)

    def _search_with_adapter_defaults(self, db_conn: sqlite3.Connection, query: str, search_type: str) -> pd.DataFrame:
        """Execute search using adapter default implementations."""
        if not self.validate_search_type(search_type):
            raise ValueError(f"Unsupported search type: {search_type}")
        
        # Route to appropriate search method
        if search_type == 'igdb_id':
            return self._search_by_igdb_id(db_conn, query)
        elif search_type == 'ens_id':
            return self._search_by_ens_id(db_conn, query)
        elif search_type == 'gene_name':
            return self._search_by_gene_name(db_conn, query)
        elif search_type == 'go_id':
            return self._search_by_go_id(db_conn, query)
        elif search_type == 'igdb_tf':
            return self._search_by_igdb_tf(db_conn, query)
        elif search_type == 'interpro_id':
            return self._search_by_interpro_id(db_conn, query)
        elif search_type.endswith('_ortho_id'):
            species = search_type.replace('_ortho_id', '').title()
            return self._search_by_ortholog_id(db_conn, query, species)
        elif search_type.endswith('_ortho_name'):
            species = search_type.replace('_ortho_name', '').title()
            return self._search_by_ortholog_name(db_conn, query, species)
        else:
            raise ValueError(f"Search type not implemented: {search_type}")
    
    def _search_by_igdb_id(self, db_conn: sqlite3.Connection, query: str) -> pd.DataFrame:
        """Search by IGDB ID and return complete ENS+IGDB information."""
        sql = """
        SELECT DISTINCT 
            l.IGDB, l.chromosome, l.start, l.end, l.strand,
            m.ENS, e.GeneName, e.GeneType, e.Description
        FROM IGDB_Locus l
        LEFT JOIN IGDB_ENS m ON l.IGDB = m.IGDB
        LEFT JOIN ENS_INFO e ON m.ENS = e.ENS
        WHERE l.IGDB = ?
        ORDER BY l.IGDB
        """
        result = pd.read_sql_query(sql, db_conn, params=[query])
        
        # If no ENS mapping found, create virtual ENS ID for medaka-specific genes
        if not result.empty and result['ENS'].isna().all():
            query_type = self.detect_id_format(query)
            virtual_ens = f"MEDAKA_{query.replace('IGDB', '')}" if query_type == 'igdb_id' else f"MEDAKA_{query}"
            result.loc[result['ENS'].isna(), 'ENS'] = virtual_ens
            
            # If no gene name, use IGDB as gene name
            if result['GeneName'].isna().all():
                result.loc[result['GeneName'].isna(), 'GeneName'] = query
                result.loc[result['GeneType'].isna(), 'GeneType'] = 'medaka_specific'
                result.loc[result['Description'].isna(), 'Description'] = f'Medaka-specific gene {query}'
        
        return result
    
    def _search_by_ens_id(self, db_conn: sqlite3.Connection, query: str) -> pd.DataFrame:
        """Search by ENS ID and return complete IGDB+ENS information (IGDB as primary)."""
        sql = """
        SELECT DISTINCT
            m.IGDB, l.chromosome, l.start, l.end, l.strand,
            e.ENS, e.GeneName, e.GeneType, e.Description
        FROM ENS_INFO e
        JOIN IGDB_ENS m ON e.ENS = m.ENS
        JOIN IGDB_Locus l ON m.IGDB = l.IGDB
        WHERE e.ENS = ?
        ORDER BY m.IGDB
        """
        return pd.read_sql_query(sql, db_conn, params=[query])
    
    def _search_by_gene_name(self, db_conn: sqlite3.Connection, query: str) -> pd.DataFrame:
        """Search by gene name in both ENS_INFO and IGDB_INFO tables (IGDB as primary)."""
        
        # First search in ENS_INFO table (original ENS gene names)
        ens_sql = """
        SELECT DISTINCT
            m.IGDB, l.chromosome, l.start, l.end, l.strand,
            e.ENS, e.GeneName, e.GeneType, e.Description
        FROM ENS_INFO e
        JOIN IGDB_ENS m ON e.ENS = m.ENS
        JOIN IGDB_Locus l ON m.IGDB = l.IGDB
        WHERE e.GeneName LIKE ? OR e.GeneName = ?
        ORDER BY m.IGDB
        """
        ens_results = pd.read_sql_query(ens_sql, db_conn, params=[f'%{query}%', query])
        
        # Then search in IGDB_NAME table (corrected IGDB gene names)
        igdb_sql = """
        SELECT DISTINCT
            l.IGDB, l.chromosome, l.start, l.end, l.strand,
            m.ENS, e.GeneName, e.GeneType, e.Description
        FROM IGDB_NAME i
        JOIN IGDB_Locus l ON i.IGDB = l.IGDB
        LEFT JOIN IGDB_ENS m ON l.IGDB = m.IGDB
        LEFT JOIN ENS_INFO e ON m.ENS = e.ENS
        WHERE i.GeneName LIKE ? OR i.GeneName = ?
        ORDER BY l.IGDB
        """
        igdb_results = pd.read_sql_query(igdb_sql, db_conn, params=[f'%{query}%', query])
        
        # Combine results and remove complete duplicates (keep multiple ENS per IGDB)
        if not ens_results.empty and not igdb_results.empty:
            combined_results = pd.concat([ens_results, igdb_results], ignore_index=True)
            return combined_results.drop_duplicates().reset_index(drop=True)
        elif not ens_results.empty:
            return ens_results
        elif not igdb_results.empty:
            return igdb_results
        else:
            return pd.DataFrame()
    
    def _search_by_go_id(self, db_conn: sqlite3.Connection, query: str) -> pd.DataFrame:
        """Search by GO ID and return complete IGDB+ENS information (IGDB as primary)."""
        sql = """
        SELECT DISTINCT
            m.IGDB, l.chromosome, l.start, l.end, l.strand,
            e.ENS, e.GeneName, e.GeneType, e.Description
        FROM ENS_INFO e
        JOIN ENS_GO g ON e.ENS = g.ENS
        JOIN IGDB_ENS m ON e.ENS = m.ENS
        JOIN IGDB_Locus l ON m.IGDB = l.IGDB
        WHERE g.GO_Accession = ?
        ORDER BY m.IGDB
        """
        return pd.read_sql_query(sql, db_conn, params=[query])
    
    def _search_by_igdb_tf(self, db_conn: sqlite3.Connection, query: str) -> pd.DataFrame:
        """Search by IGDB TF family and return complete information."""
        sql = """
        SELECT DISTINCT
            t.IGDB, l.chromosome, l.start, l.end, l.strand,
            m.ENS, e.GeneName, e.GeneType, e.Description,
            t.TF_Family
        FROM IGDB_TF t
        JOIN IGDB_Locus l ON t.IGDB = l.IGDB
        LEFT JOIN IGDB_ENS m ON t.IGDB = m.IGDB
        LEFT JOIN ENS_INFO e ON m.ENS = e.ENS
        WHERE t.TF_Family = ?
        ORDER BY t.IGDB
        """
        return pd.read_sql_query(sql, db_conn, params=[query])
    
    def _search_by_interpro_id(self, db_conn: sqlite3.Connection, query: str) -> pd.DataFrame:
        """Search by InterPro ID and return complete IGDB+ENS information (IGDB as primary)."""
        sql = """
        SELECT DISTINCT
            m.IGDB, l.chromosome, l.start, l.end, l.strand,
            e.ENS, e.GeneName, e.GeneType, e.Description
        FROM ENS_INFO e
        JOIN ENS_INTERPRO i ON e.ENS = i.ENS
        JOIN IGDB_ENS m ON e.ENS = m.ENS
        JOIN IGDB_Locus l ON m.IGDB = l.IGDB
        WHERE i.IPRID = ?
        ORDER BY m.IGDB
        """
        return pd.read_sql_query(sql, db_conn, params=[query])
    
    def _search_by_ortholog_id(self, db_conn: sqlite3.Connection, query: str, species: str) -> pd.DataFrame:
        """Search by ortholog ID and return complete IGDB+ENS information (IGDB as primary)."""
        table_name = f"ENS_{species}"
        sql = f"""
        SELECT DISTINCT
            m.IGDB, l.chromosome, l.start, l.end, l.strand,
            e.ENS, e.GeneName, e.GeneType, e.Description
        FROM ENS_INFO e
        JOIN {table_name} o ON e.ENS = o.ENS
        JOIN IGDB_ENS m ON e.ENS = m.ENS
        JOIN IGDB_Locus l ON m.IGDB = l.IGDB
        WHERE o.OrthoID = ?
        ORDER BY m.IGDB
        """
        return pd.read_sql_query(sql, db_conn, params=[query])
    
    def _search_by_ortholog_name(self, db_conn: sqlite3.Connection, query: str, species: str) -> pd.DataFrame:
        """Search by ortholog name and return complete IGDB+ENS information (IGDB as primary)."""
        table_name = f"ENS_{species}"
        sql = f"""
        SELECT DISTINCT
            m.IGDB, l.chromosome, l.start, l.end, l.strand,
            e.ENS, e.GeneName, e.GeneType, e.Description
        FROM ENS_INFO e
        JOIN {table_name} o ON e.ENS = o.ENS
        JOIN IGDB_ENS m ON e.ENS = m.ENS
        JOIN IGDB_Locus l ON m.IGDB = l.IGDB
        WHERE o.OrthoName LIKE ? OR o.OrthoName = ?
        ORDER BY m.IGDB
        """
        return pd.read_sql_query(sql, db_conn, params=[f'%{query}%', query])
    
    def get_primary_key_column(self) -> str:
        """Medaka uses IGDB as primary key for expression data."""
        return 'IGDB'
    
    def get_display_columns(self) -> List[str]:
        """Return columns to display for medaka (includes both IGDB and ENS)."""
        return ['IGDB', 'ENS', 'GeneName', 'GeneType', 'Description', 'chromosome', 'start', 'end']
    
    def format_gene_label(self, row: pd.Series) -> str:
        """Format gene label for medaka using IGDB ID."""
        igdb_id = row.get('IGDB', '')
        gene_name = row.get('GeneName', '')
        
        return f"{gene_name} ({igdb_id})"
    
    def get_common_id_column(self) -> str:
        """Medaka uses IGDB as common_id in data storage."""
        return 'IGDB'

    def parse_primary_id(self, primary_id: str) -> Dict[str, str]:
        """
        Parse primary_id into IGDB (common_id) and ENS components.
        Format: 'IGDB::ENS' or just 'IGDB' (when no ENS mapping)
        """
        if '::' in primary_id:
            parts = primary_id.split('::', 1)
            return {'igdb': parts[0], 'ens': parts[1] if len(parts) > 1 else ''}
        return {'igdb': primary_id, 'ens': ''}

    def get_common_ids(self, primary_ids: List[str]) -> List[str]:
        """
        Extract common_ids (IGDB) from primary_ids for expression data queries.
        Removes duplicates since multiple IGDB::ENS can share the same IGDB.
        """
        igdb_set = set()
        for pid in primary_ids:
            parsed = self.parse_primary_id(pid)
            if parsed['igdb']:
                igdb_set.add(parsed['igdb'])
        return list(igdb_set)

    def get_gene_sequences(self, db_conn: sqlite3.Connection, gene_id: str) -> List[Dict[str, str]]:
        """Get all transcript sequences for a gene. For medaka, use IGDB_Seq table."""
        if not gene_id:
            return []
        
        try:
            cursor = db_conn.cursor()

            parsed = self.parse_primary_id(gene_id)
            if parsed.get('igdb') and '::' in gene_id:
                igdb_id = parsed['igdb']
                display_id = gene_id
            else:
                # Use configuration-driven ID detection instead of hardcoded patterns
                id_type = self.detect_id_format(gene_id)
            
                if id_type in {'ens_id', 'secondary_id'}:  # ENS format (secondary ID)
                    # Map ENS to IGDB
                    mapping_sql = "SELECT IGDB FROM IGDB_ENS WHERE ENS = ?"
                    mapping_result = cursor.execute(mapping_sql, (gene_id,)).fetchone()
                    
                    if not mapping_result:
                        self.logger.warning("No IGDB mapping found for ENS: %s", gene_id)
                        return []
                    
                    igdb_id = mapping_result['IGDB']
                    display_id = gene_id  # Keep original ENS ID for display
                elif id_type == 'igdb_id':  # IGDB format (primary ID)
                    # Already IGDB ID
                    igdb_id = gene_id
                    display_id = gene_id
                else:
                    self.logger.warning("Unknown gene ID format: %s", gene_id)
                    return []
            
            # Get sequences using IGDB ID from IGDB_Seq table
            sql = "SELECT IGDBtrpt, nseq FROM IGDB_Seq WHERE IGDB = ? ORDER BY IGDBtrpt"
            results = cursor.execute(sql, (igdb_id,)).fetchall()
            
            sequences = []
            for row in results:
                transcript_id = row['IGDBtrpt']
                sequence = row['nseq']
                
                # Format as FASTA with both IGDB and ENS info if available
                display_id_type = self.detect_id_format(display_id)
                if display_id_type == 'ens_id':
                    fasta_header = f">{transcript_id} | {display_id} | {igdb_id} | {len(sequence)} bp"
                else:
                    fasta_header = f">{transcript_id} | {igdb_id} | {len(sequence)} bp"
                fasta_sequence = self._format_sequence_lines(sequence)
                
                sequences.append({
                    'transcript_id': transcript_id,
                    'length': len(sequence),
                    'fasta': f"{fasta_header}\n{fasta_sequence}"
                })
            
            return sequences
        except Exception as e:
            self.logger.error("Error getting sequences for %s: %s", gene_id, e)
            return []
    
    def _format_sequence_lines(self, sequence: str, line_length: int = 60) -> str:
        """Format DNA sequence into lines of specified length."""
        return '\n'.join(sequence[i:i+line_length] for i in range(0, len(sequence), line_length))
    
    def _get_igdb_gene_name(self, db_conn: sqlite3.Connection, igdb_id: str) -> str:
        """Get IGDB gene name from IGDB_NAME table."""
        if not igdb_id:
            return ""
        
        try:
            cursor = db_conn.cursor()
            cursor.execute("SELECT GeneName FROM IGDB_NAME WHERE IGDB = ?", (igdb_id,))
            result = cursor.fetchone()
            return result['GeneName'] if result else ""
        except Exception as e:
            self.logger.error("Error getting IGDB gene name for %s: %s", igdb_id, e)
            return ""
    
    def _format_combined_gene_name(self, db_conn: sqlite3.Connection, ens_name: str, igdb_id: str) -> str:
        """Format gene name as 'ens name [IGDB: igdb name]' or just 'ens name' if no IGDB name."""
        if not ens_name:
            return igdb_id if igdb_id else ""
        
        # Get IGDB gene name
        igdb_name = self._get_igdb_gene_name(db_conn, igdb_id)
        
        if igdb_name:
            return f"{ens_name} [IGDB: {igdb_name}]"
        else:
            return ens_name
    
    def format_search_results(self, raw_results: pd.DataFrame, db_conn: sqlite3.Connection = None) -> List[Dict[str, Any]]:
        """Format raw search results into standardized format for medaka (IGDB primary)."""
        if raw_results.empty:
            return []

        def _format_results(conn):
            results = []
            for _, row in raw_results.iterrows():
                igdb_id = row.get('IGDB', '')
                ens_id = row.get('ENS', '')
                ens_name = row.get('GeneName', '')

                # Create primary_id: IGDB::ENS or just IGDB if no ENS
                primary_id = f"{igdb_id}::{ens_id}" if ens_id else igdb_id

                # Format gene name as "ens name [IGDB: igdb name]" or just "ens name"
                formatted_name = self._format_combined_gene_name(conn, ens_name, igdb_id)

                result = {
                    "primary_id": primary_id,  # Unique ID: IGDB::ENS for medaka
                    "secondary_id": ens_id if ens_id else None,  # ENS for Ensembl links
                    "gene_name": formatted_name,
                    "description": row.get('Description', ''),
                    "gene_type": row.get('GeneType', ''),
                }

                # Add locus information if available
                locus = self._format_locus_string(row)
                if locus:
                    result["locus"] = locus

                # Add additional fields that might be useful
                for col in ['chromosome', 'start', 'end', 'strand']:
                    if col in row and pd.notna(row[col]):
                        result[col.lower()] = row[col]

                results.append(result)
            return results

        # Use provided connection or create new one
        if db_conn is not None:
            return _format_results(db_conn)

        # Fallback: create new connection
        db_path = Path(__file__).parent.parent.parent / "data" / self.species_key / "geneinfo.db"
        try:
            with sqlite3.connect(f"file:{db_path}?mode=ro", uri=True, check_same_thread=False) as new_conn:
                new_conn.row_factory = sqlite3.Row
                return _format_results(new_conn)
        except Exception as e:
            self.logger.warning("Could not format IGDB names (fallback to raw data): %s", e)
            # Fallback if connection fails
            results = []
            for _, row in raw_results.iterrows():
                igdb_id = row.get('IGDB', '')
                ens_id = row.get('ENS', '')
                ens_name = row.get('GeneName', '')

                primary_id = f"{igdb_id}::{ens_id}" if ens_id else igdb_id
                formatted_name = ens_name if ens_name else (igdb_id if igdb_id else ens_id)

                result = {
                    "primary_id": primary_id,
                    "secondary_id": ens_id if ens_id else None,
                    "gene_name": formatted_name,
                    "description": row.get('Description', ''),
                    "gene_type": row.get('GeneType', ''),
                }
                locus = self._format_locus_string(row)
                if locus:
                    result["locus"] = locus
                for col in ['chromosome', 'start', 'end', 'strand']:
                    if col in row and pd.notna(row[col]):
                        result[col.lower()] = row[col]
                results.append(result)
            return results
    
    def get_external_links(self, result: Dict[str, Any]) -> Dict[str, str]:
        """Generate external database links for medaka."""
        links = {}
        
        # Ensembl link using secondary_id (ENS) if available
        secondary_id = result.get('secondary_id')
        if secondary_id and secondary_id.startswith('ENSORLG'):
            # For medaka, secondary_id contains ENS IDs - always generate Ensembl links
            # Use Hugo configuration to get scientific name
            ensembl_species = self._get_ensembl_species_name()
            
            # Build Ensembl URL using dynamic species name
            base_url = f"https://ensembl.org/{ensembl_species}/Gene/Summary?g="
            links["Ensembl"] = f"{base_url}{secondary_id}"
        
        return links
    
    def get_tf_families(self, db_conn: sqlite3.Connection, primary_id: str, secondary_id: str = None) -> List[str]:
        """Get transcription factor families for medaka (IGDB-based with ENS fallback)."""
        try:
            cursor = db_conn.cursor()
            
            # For medaka, try IGDB_TF table first with primary_id (IGDB)
            primary_id_type = self.detect_id_format(primary_id)
            if primary_id_type == 'igdb_id':
                cursor.execute("SELECT IGDB, TF_Family FROM IGDB_TF WHERE IGDB = ?", (primary_id,))
                igdb_tf_results = cursor.fetchall()
                if igdb_tf_results:
                    return [row['TF_Family'] for row in igdb_tf_results]
            
            # Fallback to ENS_TF table using secondary_id (ENS)
            secondary_id_type = self.detect_id_format(secondary_id) if secondary_id else 'unknown'
            if secondary_id and secondary_id_type == 'ens_id':
                cursor.execute("SELECT ENS, TF_Family FROM ENS_TF WHERE ENS = ?", (secondary_id,))
                ens_tf_results = cursor.fetchall()
                if ens_tf_results:
                    return [row['TF_Family'] for row in ens_tf_results]
                    
        except Exception as e:
            self.logger.error("Error getting TF families for medaka %s: %s", primary_id, e)
        
        return []
    
    def get_summary_columns(self) -> Dict[str, str]:
        """Return summary table column configuration for medaka."""
        return {
            'gene_name': 'Gene Name',
            'primary_id': 'IGDB ID',
            'secondary_id': 'Ensembl ID',
            'locus': 'Locus',
            'tf_families': 'TF Family',
            'description': 'Description'
        }
    
    def format_search_result_display(self, result: Dict[str, Any]) -> Dict[str, Any]:
        """Format search result display for medaka."""
        gene_name = result.get("gene_name", "N/A")
        primary_id = result.get("primary_id", "")  # IGDB::ENS or IGDB

        # Consistent with standard species: main line shows gene name, secondary line shows parsed primary_id (IGDB for medaka)
        parsed = self.parse_primary_id(primary_id) if primary_id else {"igdb": "", "ens": ""}
        common_id = parsed.get("igdb") or primary_id
        
        return {
            "main_text": gene_name,
            "secondary_items": [
                {"text": common_id, "className": "text-muted"} if common_id else None
            ]
        }
    
    def is_tf_search_field(self, field: str) -> bool:
        """Check if a search field is for transcription factor searches."""
        return field in ['tf_family', 'igdb_tf']
    
    def convert_ids_to_gene_names(self, gene_ids: List[str]) -> List[str]:
        """Convert gene IDs to gene names for medaka."""
        # Import here to avoid circular dependency
        from ..gene_search import GeneSearchAPI
        
        try:
            gene_api = GeneSearchAPI(self.species_key)
            converted = []
            
            for gene_id in gene_ids:
                # Use configuration-driven ID detection instead of hardcoded patterns
                id_type = self.detect_id_format(gene_id)
                
                if id_type in ['igdb_id', 'ens_id']:  # Recognizable ID formats
                    # Get gene details to find the gene name
                    details = gene_api.get_gene_details([gene_id])
                    if gene_id in details and 'gene_name' in details[gene_id]:
                        gene_name = details[gene_id]['gene_name']
                        if gene_name:
                            converted.append(gene_name)
                        else:
                            converted.append(gene_id)  # Keep original if no name
                    else:
                        converted.append(gene_id)  # Keep original if not found
                else:
                    converted.append(gene_id)  # Already a gene name or unknown format
            
            return converted
            
        except Exception as e:
            self.logger.error("Error converting gene IDs: %s", e)
            return gene_ids  # Return original list if conversion fails
    
    def get_column_display_config(self) -> Dict[str, str]:
        """Return column display configuration for medaka."""
        return {
            'gene_name': 'Gene Name',
            'primary_id': 'IGDB ID',
            'secondary_id': 'Ensembl ID',
            'locus': 'Locus',
            'description': 'Description',
            'gene_type': 'Gene Type',
            'tf_families': 'TF Family'
        }
