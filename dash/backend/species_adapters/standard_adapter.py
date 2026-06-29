"""
Standard species adapter for ENS-based species (zebrafish, etc.).

This adapter handles species that use the standard Ensembl data model
without additional gene ID mapping layers.
"""

import pandas as pd
import sqlite3
from typing import List, Dict, Any
from .base_adapter import SpeciesAdapter


class StandardAdapter(SpeciesAdapter):
    """Adapter for standard ENS-based species."""
    
    def get_search_fields(self) -> List[Dict[str, str]]:
        """Return search fields for standard ENS species."""
        return [
            {'label': 'Gene Name', 'value': 'gene_name'},
            {'label': 'Ensembl ID', 'value': 'ens_id'}, 
            {'label': 'GO ID', 'value': 'go_id'},
            {'label': 'TF Family', 'value': 'tf_family'},
            {'label': 'InterPro ID', 'value': 'interpro_id'},
            {'label': 'Human Ortholog ID', 'value': 'human_ortho_id'},
            {'label': 'Mouse Ortholog ID', 'value': 'mouse_ortho_id'},
            {'label': 'Human Ortholog Name', 'value': 'human_ortho_name'},
            {'label': 'Mouse Ortholog Name', 'value': 'mouse_ortho_name'},
        ]
    
    def search_genes(self, db_conn: sqlite3.Connection, query: str, search_type: str) -> pd.DataFrame:
        """
        Execute gene search using TOML config or adapter defaults.
        
        This method first checks if TOML configuration exists for the search type,
        and if so, uses that. Otherwise falls back to adapter defaults.
        """
        if not query.strip():
            return pd.DataFrame()
        
        cleaned_query = self._clean_query(query)
        if not cleaned_query:
            return pd.DataFrame()
        
        # Check if we have TOML configuration for this search type
        if search_type in self._search_field_configs:
            return self._search_with_toml_config(db_conn, cleaned_query, search_type)
        
        # Fallback to adapter defaults
        return self._search_with_adapter_defaults(db_conn, cleaned_query, search_type)
    
    def _search_with_toml_config(self, db_conn: sqlite3.Connection, query: str, search_type: str) -> pd.DataFrame:
        """Execute search using TOML configuration."""
        config = self._search_field_configs[search_type]
        db_table = config['db_table']
        db_column = config['db_column']
        match_type = config.get('match_type', 'like')
        
        # Build base query for ENS-based species
        base_select = """
        SELECT DISTINCT e.ENS, e.GeneName, e.GeneType, e.Description,
               l.chromosome, l.start, l.end, l.strand
        FROM ENS_INFO e
        LEFT JOIN ENS_Locus l ON e.ENS = l.ENS
        """
        
        # Prepare search terms
        search_terms = [term.strip() for term in query.split() if term.strip()]
        all_results = pd.DataFrame()
        
        for term in search_terms:
            operator = "=" if match_type == "exact" else "LIKE"
            search_pattern = term if match_type == "exact" else f"%{term}%"
            
            if db_table == "ENS_INFO":
                sql = f"{base_select} WHERE e.{db_column} {operator} ?"
            else:
                sql = f"{base_select} JOIN {db_table} t ON e.ENS = t.ENS WHERE t.{db_column} {operator} ?"
            
            sql += " ORDER BY e.GeneName"
            
            try:
                result_df = pd.read_sql_query(sql, db_conn, params=[search_pattern])
                if not result_df.empty:
                    all_results = pd.concat([all_results, result_df], ignore_index=True)
            except Exception as e:
                self.logger.warning("TOML search error for term '%s': %s", term, e)
        
        # Remove duplicates and return
        if not all_results.empty:
            return all_results.drop_duplicates(subset=['ENS']).reset_index(drop=True)
        return pd.DataFrame()
    
    def _search_with_adapter_defaults(self, db_conn: sqlite3.Connection, query: str, search_type: str) -> pd.DataFrame:
        """Execute search using adapter default implementations."""
        if not self.validate_search_type(search_type):
            raise ValueError(f"Unsupported search type: {search_type}")
        
        # Route to appropriate search method
        if search_type == 'gene_name':
            return self._search_by_gene_name(db_conn, query)
        elif search_type == 'ens_id':
            return self._search_by_ens_id(db_conn, query)
        elif search_type == 'go_id':
            return self._search_by_go_id(db_conn, query)
        elif search_type == 'tf_family':
            return self._search_by_tf_family(db_conn, query)
        elif search_type == 'interpro_id':
            return self._search_by_interpro_id(db_conn, query)
        elif search_type == 'human_ortho_id':
            return self._search_by_ortholog_id(db_conn, query, 'Human')
        elif search_type == 'mouse_ortho_id':
            return self._search_by_ortholog_id(db_conn, query, 'Mouse')
        elif search_type == 'human_ortho_name':
            return self._search_by_ortholog_name(db_conn, query, 'Human')
        elif search_type == 'mouse_ortho_name':
            return self._search_by_ortholog_name(db_conn, query, 'Mouse')
        else:
            raise ValueError(f"Search type not implemented: {search_type}")
    
    def batch_search_genes(self, db_conn: sqlite3.Connection, gene_ids: List[str], search_type: str) -> pd.DataFrame:
        """
        Execute gene search for multiple genes efficiently using IN clause for standard species.
        
        This method is optimized for 'ens_id' or 'primary_id' search types.
        For other types, it falls back to individual searches or raises an error.
        """
        if not gene_ids:
            return pd.DataFrame()

        placeholders = ','.join(['?' for _ in gene_ids])
        
        if search_type == 'ens_id' or search_type == 'primary_id':
            sql = f"""
            SELECT DISTINCT e.ENS, e.GeneName, e.GeneType, e.Description,
                   l.chromosome, l.start, l.end, l.strand
            FROM ENS_INFO e
            LEFT JOIN ENS_Locus l ON e.ENS = l.ENS
            WHERE e.ENS IN ({placeholders})
            ORDER BY e.ENS
            """
            return pd.read_sql_query(sql, db_conn, params=gene_ids)
        elif search_type == 'gene_name':
            # Batch search by name is complex due to LIKE operator.
            # For now, if exact name match is needed, run individually.
            # If partial match, a single query with multiple LIKEs can be built
            # but usually for N+1 it's about exact ID lookups.
            # Fallback to individual search for now to avoid complexity or incorrect behavior.
            self.logger.warning(f"Batch search for 'gene_name' is not optimized; performing individual searches. Consider using ID-based search for performance.")
            all_results = pd.DataFrame()
            for gene_id in gene_ids: # Note: gene_ids here are treated as individual gene names
                result = self._search_by_gene_name(db_conn, gene_id)
                if not result.empty:
                    all_results = pd.concat([all_results, result], ignore_index=True)
            return all_results.drop_duplicates(subset=['ENS']).reset_index(drop=True)
        else:
            # For other fields, a generic batch query might be possible but would require
            # building dynamic SQL based on search_type, which is out of scope for simple N+1 fix.
            self.logger.warning(f"Batch search type '{search_type}' not directly optimized for standard species; performing individual searches.")
            all_results = pd.DataFrame()
            for gene_id in gene_ids:
                result = self.search_genes(db_conn, gene_id, search_type) # Calls the single-query search_genes
                if not result.empty:
                    all_results = pd.concat([all_results, result], ignore_index=True)
            return all_results.drop_duplicates(subset=['ENS']).reset_index(drop=True)

    def _search_by_gene_name(self, db_conn: sqlite3.Connection, query: str) -> pd.DataFrame:
        """Search by gene name in ENS_INFO table."""
        sql = """
        SELECT DISTINCT e.ENS, e.GeneName, e.GeneType, e.Description,
               l.chromosome, l.start, l.end, l.strand
        FROM ENS_INFO e
        LEFT JOIN ENS_Locus l ON e.ENS = l.ENS
        WHERE e.GeneName LIKE ? OR e.GeneName = ?
        ORDER BY e.ENS
        """
        return pd.read_sql_query(sql, db_conn, params=[f'%{query}%', query])
    
    def _search_by_ens_id(self, db_conn: sqlite3.Connection, query: str) -> pd.DataFrame:
        """Search by Ensembl ID."""
        sql = """
        SELECT DISTINCT e.ENS, e.GeneName, e.GeneType, e.Description,
               l.chromosome, l.start, l.end, l.strand
        FROM ENS_INFO e
        LEFT JOIN ENS_Locus l ON e.ENS = l.ENS
        WHERE e.ENS = ?
        ORDER BY e.ENS
        """
        return pd.read_sql_query(sql, db_conn, params=[query])
    
    def _search_by_go_id(self, db_conn: sqlite3.Connection, query: str) -> pd.DataFrame:
        """Search by GO ID and return corresponding ENS genes."""
        sql = """
        SELECT DISTINCT e.ENS, e.GeneName, e.GeneType, e.Description,
               l.chromosome, l.start, l.end, l.strand
        FROM ENS_INFO e
        JOIN ENS_GO g ON e.ENS = g.ENS
        LEFT JOIN ENS_Locus l ON e.ENS = l.ENS
        WHERE g.GO_Accession = ?
        ORDER BY e.ENS
        """
        return pd.read_sql_query(sql, db_conn, params=[query])
    
    def _search_by_tf_family(self, db_conn: sqlite3.Connection, query: str) -> pd.DataFrame:
        """Search by transcription factor family."""
        sql = """
        SELECT DISTINCT e.ENS, e.GeneName, e.GeneType, e.Description,
               l.chromosome, l.start, l.end, l.strand
        FROM ENS_INFO e
        JOIN ENS_TF t ON e.ENS = t.ENS
        LEFT JOIN ENS_Locus l ON e.ENS = l.ENS
        WHERE t.TF_Family = ?
        ORDER BY e.ENS
        """
        return pd.read_sql_query(sql, db_conn, params=[query])
    
    def _search_by_interpro_id(self, db_conn: sqlite3.Connection, query: str) -> pd.DataFrame:
        """Search by InterPro ID."""
        sql = """
        SELECT DISTINCT e.ENS, e.GeneName, e.GeneType, e.Description,
               l.chromosome, l.start, l.end, l.strand
        FROM ENS_INFO e
        JOIN ENS_INTERPRO i ON e.ENS = i.ENS
        LEFT JOIN ENS_Locus l ON e.ENS = l.ENS
        WHERE i.IPRID = ?
        ORDER BY e.ENS
        """
        return pd.read_sql_query(sql, db_conn, params=[query])
    
    def _search_by_ortholog_id(self, db_conn: sqlite3.Connection, query: str, species: str) -> pd.DataFrame:
        """Search by ortholog ID from specified species."""
        table_name = f"ENS_Ortho{species}"
        sql = f"""
        SELECT DISTINCT e.ENS, e.GeneName, e.GeneType, e.Description,
               l.chromosome, l.start, l.end, l.strand
        FROM ENS_INFO e
        JOIN {table_name} o ON e.ENS = o.ENS
        LEFT JOIN ENS_Locus l ON e.ENS = l.ENS
        WHERE o.OrthoID = ?
        ORDER BY e.ENS
        """
        return pd.read_sql_query(sql, db_conn, params=[query])
    
    def _search_by_ortholog_name(self, db_conn: sqlite3.Connection, query: str, species: str) -> pd.DataFrame:
        """Search by ortholog name from specified species."""
        table_name = f"ENS_Ortho{species}"
        sql = f"""
        SELECT DISTINCT e.ENS, e.GeneName, e.GeneType, e.Description,
               l.chromosome, l.start, l.end, l.strand
        FROM ENS_INFO e
        JOIN {table_name} o ON e.ENS = o.ENS
        LEFT JOIN ENS_Locus l ON e.ENS = l.ENS
        WHERE o.OrthoName LIKE ? OR o.OrthoName = ?
        ORDER BY e.ENS
        """
        return pd.read_sql_query(sql, db_conn, params=[f'%{query}%', query])
    
    def get_primary_key_column(self) -> str:
        """Standard species use ENS as primary key."""
        return 'ENS'
    
    def get_display_columns(self) -> List[str]:
        """Return columns to display for standard species."""
        return ['ENS', 'GeneName', 'GeneType', 'Description']
    
    def format_gene_label(self, row: pd.Series) -> str:
        """Format gene label for standard species."""
        ens_id = row.get('ENS', '')
        gene_name = row.get('GeneName', '')
        
        # Use last 5 characters of ENS ID for compact display
        if len(ens_id) > 5:
            ens_short = ens_id[-5:]
        else:
            ens_short = ens_id
            
        return f"{gene_name} ({ens_short})"
    
    def get_common_id_column(self) -> str:
        """Standard species use ENS as common_id in data storage."""
        return 'ENS'

    def parse_primary_id(self, primary_id: str) -> Dict[str, str]:
        """For standard species, primary_id equals ENS ID (common_id)."""
        return {'ens': primary_id, 'igdb': ''}

    def get_common_ids(self, primary_ids: List[str]) -> List[str]:
        """For standard species, primary_ids are already common_ids (ENS)."""
        return [pid for pid in primary_ids if pid]

    def get_gene_sequences(self, db_conn, gene_id: str) -> List[Dict[str, str]]:
        """Get all transcript sequences for a gene from ENS_Seq table."""
        if not gene_id:
            return []
        
        try:
            cursor = db_conn.cursor()
            sql = "SELECT ENStrpt, nseq FROM ENS_Seq WHERE ENS = ? ORDER BY ENStrpt"
            results = cursor.execute(sql, (gene_id,)).fetchall()
            
            sequences = []
            for row in results:
                transcript_id = row['ENStrpt']
                sequence = row['nseq']
                
                # Format as FASTA
                fasta_header = f">{transcript_id} | {gene_id} | {len(sequence)} bp"
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
    
    def format_search_results(self, raw_results: pd.DataFrame, db_conn=None) -> List[Dict[str, Any]]:
        """Format raw search results into standardized format for standard species."""
        if raw_results.empty:
            return []

        results = []
        for _, row in raw_results.iterrows():
            primary_id = row.get('ENS', '')
            result = {
                "primary_id": primary_id,  # ENS is primary for standard species
                "secondary_id": None,  # No secondary ID for standard species
                "gene_name": row.get('GeneName', ''),
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
    
    def get_external_links(self, result: Dict[str, Any]) -> Dict[str, str]:
        """Generate external database links for standard species."""
        links = {}
        
        # Ensembl link for primary ID if it's an ENS format
        primary_id = result.get('primary_id')
        if primary_id:
            id_type = self.detect_id_format(primary_id)
            if id_type == 'ens_id':
                # Use Hugo configuration to get scientific name
                ensembl_species = self._get_ensembl_species_name()
                
                # Build Ensembl URL using dynamic species name
                base_url = f"https://ensembl.org/{ensembl_species}/Gene/Summary?g="
                links["Ensembl"] = f"{base_url}{primary_id}"
        
        return links
    
    def get_tf_families(self, db_conn, primary_id: str, secondary_id: str = None) -> List[str]:
        """Get transcription factor families for standard species (ENS-based)."""
        try:
            cursor = db_conn.cursor()
            # For standard species, primary_id is ENS
            cursor.execute("SELECT ENS, TF_Family FROM ENS_TF WHERE ENS = ?", (primary_id,))
            tf_results = cursor.fetchall()
            if tf_results:
                return [row['TF_Family'] for row in tf_results]
        except Exception as e:
            self.logger.error("Error getting TF families for %s: %s", primary_id, e)
        
        return []
    
    def get_summary_columns(self) -> Dict[str, str]:
        """Return summary table column configuration for standard species."""
        return {
            'gene_name': 'Gene Name',
            'primary_id': 'Ensembl ID',
            'locus': 'Locus',
            'tf_families': 'TF Family',
            'description': 'Description'
        }
    
    def format_search_result_display(self, result: Dict[str, Any]) -> Dict[str, Any]:
        """Format search result display for standard species."""
        gene_name = result.get("gene_name", "N/A")
        primary_id = result.get("primary_id", "")
        
        return {
            "main_text": gene_name,
            "secondary_items": [
                {"text": primary_id, "className": "text-muted"} if primary_id else None
            ]
        }
    
    def is_tf_search_field(self, field: str) -> bool:
        """Check if a search field is for transcription factor searches."""
        return field in ['tf_family']
    
    def convert_ids_to_gene_names(self, gene_ids: List[str]) -> List[str]:
        """Convert gene IDs to gene names for standard species."""
        # Import here to avoid circular dependency
        from ..gene_search import GeneSearchAPI
        
        try:
            gene_api = GeneSearchAPI(self.species_key)
            converted = []
            
            for gene_id in gene_ids:
                # Use configuration-driven ID detection instead of hardcoded patterns
                id_type = self.detect_id_format(gene_id)
                
                if id_type in ['ens_id', 'primary_id']:  # Recognizable ID formats
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
        """Return column display configuration for standard species."""
        return {
            'gene_name': 'Gene Name',
            'primary_id': 'Ensembl ID',
            'secondary_id': 'Secondary ID',  # Unused for standard species
            'locus': 'Locus',
            'description': 'Description',
            'gene_type': 'Gene Type',
            'tf_families': 'TF Family'
        }
