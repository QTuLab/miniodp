"""
Gene Search API - Pure Adapter Pattern Implementation

This module provides intelligent gene search functionality through a clean
adapter pattern that eliminates all species-specific code from the API layer.
"""

import sqlite3
import contextlib
import logging
from typing import Any, Dict, Iterator, List, Optional, Sequence
from pathlib import Path
from .species_adapters import create_species_adapter


logger = logging.getLogger(__name__)


class GeneSearchAPI:
    """Pure adapter-based gene search API with zero species-specific code."""
    
    def __init__(self, species_key: str):
        if not species_key:
            raise ValueError("Species key is required and cannot be empty")
        if not isinstance(species_key, str):
            raise TypeError("Species key must be a string")
        
        self.species = species_key.strip()
        if not self.species:
            raise ValueError("Species key cannot be empty or whitespace")
            
        # Use global DASH_DATA_PATH or fall back to legacy path
        import os
        data_root = os.environ.get('DASH_DATA_PATH', Path(__file__).parent.parent / "data")
        self.data_dir = Path(data_root) / self.species
        self.db_path = self.data_dir / "geneinfo.db"
        
        # Initialize species adapter (handles all species-specific logic)
        self.adapter = create_species_adapter(self.species)
        
        if not self.db_path.exists():
            raise FileNotFoundError(f"Gene info database not found: {self.db_path}")

    @contextlib.contextmanager
    def _get_connection(self) -> Iterator[sqlite3.Connection]:
        """Get database connection with read-only access."""
        conn = sqlite3.connect(f"file:{self.db_path}?mode=ro", uri=True, check_same_thread=False)
        try:
            conn.row_factory = sqlite3.Row
            conn.execute("PRAGMA query_only = ON")
            yield conn
        finally:
            conn.close()

    def _normalize_gene_ids(self, gene_ids: Optional[Sequence[str]]) -> List[str]:
        """Ensure provided gene IDs are non-empty strings."""
        if gene_ids is None:
            return []
        if isinstance(gene_ids, str):
            raise TypeError("gene_ids must be a sequence of strings, not a single string")
        normalized: List[str] = []
        for gene_id in gene_ids:
            if gene_id is None:
                continue
            if not isinstance(gene_id, str):
                raise TypeError(f"Gene ID must be a string, got {type(gene_id)!r}")
            cleaned = gene_id.strip()
            if cleaned:
                normalized.append(cleaned)
        return normalized

    def search_genes(self, query: str, field: str = 'gene_name', limit: int = 50) -> List[Dict[str, Any]]:
        """
        Search genes using pure adapter pattern.
        
        This method delegates all search logic to the species adapter,
        ensuring zero species-specific code in the API layer.
        
        Args:
            query: Search query string
            field: Search field type
            limit: Maximum number of results
            
        Returns:
            List of standardized search results
        """
        if query is None:
            return []
        if not isinstance(query, str):
            raise TypeError("query must be a string")
        query = query.strip()
        if not query:
            return []

        if field is None:
            field = 'gene_name'
        if not isinstance(field, str):
            raise TypeError("field must be a string")
        field = field.strip() or 'gene_name'

        if not isinstance(limit, int):
            raise TypeError("limit must be an integer")
        if limit <= 0:
            raise ValueError("limit must be greater than zero")
        
        # Validate search field
        if not self.adapter.validate_search_type(field):
            logger.warning("Unsupported search type '%s' for %s", field, self.species)
            return []
        
        try:
            with self._get_connection() as conn:
                # Adapter handles all search logic (TOML config + species logic)
                raw_results = self.adapter.search_genes(conn, query, field)

                # Adapter formats results into standardized format (pass conn for efficiency)
                formatted_results = self.adapter.format_search_results(raw_results, db_conn=conn)

                # Add backward compatibility fields for UI
                for result in formatted_results:
                    self._add_backward_compatibility_fields(result)

                # Apply limit
                return formatted_results[:limit]

        except Exception as e:
            logger.error("Search error for %s: %s", self.species, e)
            return []

    def get_gene_names_by_ids(self, gene_ids: Sequence[str]) -> Dict[str, str]:
        """
        Get gene names by gene IDs using adapter pattern, optimized for batch queries.
        Supports composite IDs (e.g., 'IGDB::ENS' for medaka).

        Args:
            gene_ids: List of gene identifiers (can be composite IDs)

        Returns:
            Dict mapping gene IDs to gene names
        """
        normalized_ids = self._normalize_gene_ids(gene_ids)
        if not normalized_ids:
            return {}

        result_map = {gid: gid for gid in normalized_ids}  # Initialize with self-mapping as fallback

        # Separate primary_ids with :: from simple IDs
        # common_id_map: {common_id(IGDB): [(primary_id, ens), ...]}
        common_id_map: Dict[str, List[tuple]] = {}
        simple_ids: Dict[str, List[str]] = {}

        for gene_id in normalized_ids:
            parsed = self.adapter.parse_primary_id(gene_id)
            if parsed.get('igdb') and '::' in gene_id:
                common_id = parsed['igdb']
                ens = parsed.get('ens', '')
                if common_id not in common_id_map:
                    common_id_map[common_id] = []
                common_id_map[common_id].append((gene_id, ens))
            else:
                detected_type = self._detect_id_type(gene_id)
                if detected_type not in simple_ids:
                    simple_ids[detected_type] = []
                simple_ids[detected_type].append(gene_id)

        try:
            with self._get_connection() as conn:
                # Handle primary_ids with :: (medaka IGDB::ENS format)
                if common_id_map:
                    common_id_list = list(common_id_map.keys())
                    raw_results = self.adapter.batch_search_genes(conn, common_id_list, 'igdb_id')
                    formatted_results = self.adapter.format_search_results(raw_results, db_conn=conn)

                    for res in formatted_results:
                        # Parse primary_id to get common_id (IGDB)
                        res_parsed = self.adapter.parse_primary_id(res.get('primary_id', ''))
                        common_id = res_parsed.get('igdb', '')
                        ens = res.get('secondary_id', '') or ''
                        gene_name = res.get('gene_name', common_id)

                        if common_id in common_id_map:
                            for primary_id, expected_ens in common_id_map[common_id]:
                                if ens == expected_ens:
                                    result_map[primary_id] = gene_name

                # Handle simple IDs (original logic)
                for search_type, ids_to_batch in simple_ids.items():
                    if search_type in ['ens_id', 'igdb_id', 'primary_id', 'secondary_id']:
                        raw_results = self.adapter.batch_search_genes(conn, ids_to_batch, search_type)
                        formatted_results = self.adapter.format_search_results(raw_results, db_conn=conn)

                        for res in formatted_results:
                            original_queried_id = None
                            if search_type == 'igdb_id' and res.get('primary_id') in ids_to_batch:
                                original_queried_id = res['primary_id']
                            elif search_type == 'ens_id' and res.get('secondary_id') in ids_to_batch:
                                original_queried_id = res['secondary_id']
                            elif search_type == 'ens_id' and res.get('primary_id') in ids_to_batch:
                                original_queried_id = res['primary_id']

                            if original_queried_id:
                                result_map[original_queried_id] = res.get('gene_name') or original_queried_id

                    else:
                        for gene_id in ids_to_batch:
                            raw_results = self.adapter.search_genes(conn, gene_id, search_type)
                            formatted_results = self.adapter.format_search_results(raw_results, db_conn=conn)
                            if formatted_results:
                                result_map[gene_id] = formatted_results[0].get('gene_name') or gene_id
        except Exception as e:
            logger.error("Error getting gene names for %s: %s", self.species, e)

        return result_map

    def get_gene_details(self, gene_ids: Sequence[str]) -> Dict[str, Dict[str, Any]]:
        """
        Get detailed gene information using adapter pattern, optimized for batch queries.
        Supports primary_ids (e.g., 'IGDB::ENS' for medaka) for unique identification.

        Args:
            gene_ids: List of gene identifiers (can be primary_ids with :: format)

        Returns:
            Dict mapping gene IDs to detailed information
        """
        normalized_ids = self._normalize_gene_ids(gene_ids)
        if not normalized_ids:
            return {}

        details: Dict[str, Dict[str, Any]] = {gid: {"gene_id": gid, "gene_name": gid} for gid in normalized_ids}

        # Parse primary_ids and group by common_id (IGDB) for batch search
        # common_id_map: {common_id: [(primary_id, ens_id), ...]}
        common_id_map: Dict[str, List[tuple]] = {}
        simple_ids: Dict[str, List[str]] = {}  # Non-composite IDs grouped by type

        for gene_id in normalized_ids:
            parsed = self.adapter.parse_primary_id(gene_id)
            if parsed.get('igdb') and '::' in gene_id:
                # This is a primary_id with IGDB::ENS format
                common_id = parsed['igdb']
                ens = parsed.get('ens', '')
                if common_id not in common_id_map:
                    common_id_map[common_id] = []
                common_id_map[common_id].append((gene_id, ens))
            else:
                # Simple ID - use existing logic
                detected_type = self._detect_id_type(gene_id)
                if detected_type not in simple_ids:
                    simple_ids[detected_type] = []
                simple_ids[detected_type].append(gene_id)

        try:
            with self._get_connection() as conn:
                cursor = conn.cursor()

                # Handle primary_ids with :: (search by common_id/IGDB, filter by ENS)
                if common_id_map:
                    common_id_list = list(common_id_map.keys())
                    raw_results = self.adapter.batch_search_genes(conn, common_id_list, 'igdb_id')
                    formatted_results = self.adapter.format_search_results(raw_results, db_conn=conn)

                    for gene_info in formatted_results:
                        # Parse primary_id to get common_id (IGDB)
                        res_parsed = self.adapter.parse_primary_id(gene_info.get('primary_id', ''))
                        common_id = res_parsed.get('igdb', '')
                        ens = gene_info.get('secondary_id', '') or ''

                        # Find matching primary_id
                        if common_id in common_id_map:
                            for primary_id, expected_ens in common_id_map[common_id]:
                                if ens == expected_ens:
                                    self._add_backward_compatibility_fields(gene_info)
                                    gene_info['external_links'] = self.adapter.get_external_links(gene_info)
                                    self._add_comprehensive_details(cursor, gene_info, common_id)
                                    details[primary_id] = gene_info.copy()
                                    break

                # Handle simple IDs (original logic)
                for search_type, ids_to_batch in simple_ids.items():
                    if search_type in ['ens_id', 'igdb_id', 'primary_id', 'secondary_id']:
                        raw_results = self.adapter.batch_search_genes(conn, ids_to_batch, search_type)
                        formatted_results = self.adapter.format_search_results(raw_results, db_conn=conn)

                        for gene_info in formatted_results:
                            original_queried_id = None
                            if search_type == 'igdb_id' and gene_info.get('primary_id') in ids_to_batch:
                                original_queried_id = gene_info['primary_id']
                            elif search_type == 'ens_id' and gene_info.get('secondary_id') in ids_to_batch:
                                original_queried_id = gene_info['secondary_id']
                            elif search_type == 'ens_id' and gene_info.get('primary_id') in ids_to_batch:
                                original_queried_id = gene_info['primary_id']

                            if original_queried_id:
                                self._add_backward_compatibility_fields(gene_info)
                                gene_info['external_links'] = self.adapter.get_external_links(gene_info)
                                self._add_comprehensive_details(cursor, gene_info, gene_info.get('primary_id', original_queried_id))
                                details[original_queried_id] = gene_info

                    else:
                        for gene_id in ids_to_batch:
                            raw_results = self.adapter.search_genes(conn, gene_id, search_type)
                            formatted_results = self.adapter.format_search_results(raw_results, db_conn=conn)

                            if formatted_results:
                                gene_info = formatted_results[0]
                                self._add_backward_compatibility_fields(gene_info)
                                gene_info['external_links'] = self.adapter.get_external_links(gene_info)
                                self._add_comprehensive_details(cursor, gene_info, gene_info.get('primary_id', gene_id))
                                details[gene_id] = gene_info

        except Exception as e:
            logger.error("Error getting gene details for %s: %s", self.species, e)

        return details

    def get_gene_sequence(self, gene_id: str) -> List[Dict[str, str]]:
        """
        Get all transcript sequences for a gene using adapter pattern.
        
        Args:
            gene_id: Gene identifier
            
        Returns:
            List of transcript sequence dictionaries
        """
        if not gene_id:
            return []
        
        try:
            with self._get_connection() as conn:
                return self.adapter.get_gene_sequences(conn, gene_id)
        except Exception as e:
            logger.error("Error getting gene sequence for %s: %s", self.species, e)
            return []

    def get_available_search_fields(self) -> List[Dict[str, str]]:
        """
        Get available search fields from adapter.
        
        Returns:
            List of search field configurations
        """
        return self.adapter.get_available_search_fields()

    def _detect_id_type(self, gene_id: str) -> str:
        """
        Detect the type of gene ID and return appropriate search field.
        Uses adapter's configuration-driven detection method.
        
        Args:
            gene_id: Gene identifier
            
        Returns:
            Search field type or 'gene_name' if unknown
        """
        if not isinstance(gene_id, str):
            raise TypeError("gene_id must be a string")
        cleaned = gene_id.strip()
        if not cleaned:
            return 'gene_name'
        return self.adapter.detect_id_format(cleaned)

    def _add_comprehensive_details(self, cursor: sqlite3.Cursor, gene_info: Dict[str, Any], primary_id: str) -> None:
        """
        Add comprehensive gene details using database queries.
        
        This method adds GO terms, InterPro, orthologs, and TF families
        to the gene information dictionary.
        """
        # Use adapter to get the appropriate ID for ENS-based queries
        # For standard species: primary_id = ENS, secondary_id = None
        # For medaka: primary_id = IGDB, secondary_id = ENS
        ens_id = gene_info.get('secondary_id') if gene_info.get('secondary_id') else primary_id
        
        # GO Terms
        try:
            cursor.execute(
                "SELECT ENS, GO_Accession, GO_Name, GO_Domain FROM ENS_GO WHERE ENS = ? ORDER BY GO_Domain, GO_Accession",
                (ens_id,)
            )
            go_results = cursor.fetchall()
            if go_results:
                gene_info['go_terms'] = [dict(row) for row in go_results]
        except Exception as go_error:
            logger.debug("GO lookup skipped for %s: %s", ens_id, go_error)

        # InterPro
        try:
            cursor.execute(
                "SELECT ENS, IPRID, IPRShortDesc FROM ENS_INTERPRO WHERE ENS = ? ORDER BY IPRID",
                (ens_id,)
            )
            ipr_results = cursor.fetchall()
            if ipr_results:
                gene_info['interpro'] = [dict(row) for row in ipr_results]
        except Exception as interpro_error:
            logger.debug("InterPro lookup skipped for %s: %s", ens_id, interpro_error)

        # Orthologs
        for ortho_species in ['Human', 'Mouse']:
            try:
                cursor.execute(
                    f"SELECT ENS, OrthoID, OrthoName FROM ENS_Ortho{ortho_species} WHERE ENS = ?",
                    (ens_id,)
                )
                ortho_results = cursor.fetchall()
                if ortho_results:
                    key = f"orthologs_{ortho_species.lower()}"
                    gene_info[key] = [dict(row) for row in ortho_results]
            except Exception as ortho_error:
                logger.debug(
                    "Ortholog lookup for %s (%s) failed: %s",
                    ens_id,
                    ortho_species,
                    ortho_error,
                )

        # Transcription Factor Information - pure adapter delegation
        try:
            with self._get_connection() as tf_conn:
                tf_families = self.adapter.get_tf_families(
                    tf_conn, 
                    primary_id, 
                    gene_info.get('secondary_id')
                )
                if tf_families:
                    gene_info['tf_families'] = tf_families
        except Exception as tf_error:
            logger.debug(
                "Transcription factor lookup skipped for %s: %s", primary_id, tf_error
            )

    def _add_backward_compatibility_fields(self, gene_info: Dict[str, Any]) -> None:
        """
        Add backward compatibility fields for UI components.
        
        Maps the new standardized format (primary_id/secondary_id) to 
        the legacy format (igdb_id/ens_id) that the UI expects.
        Uses adapter configuration to determine ID types.
        """
        primary_id = gene_info.get('primary_id', '')
        secondary_id = gene_info.get('secondary_id')
        
        # Use adapter to detect ID types instead of hardcoded patterns
        primary_type = self.adapter.detect_id_format(primary_id) if primary_id else 'unknown'
        secondary_type = self.adapter.detect_id_format(secondary_id) if secondary_id else 'unknown'
        
        # Map to legacy format based on detected types
        if primary_type == 'igdb_id':
            # Species with IGDB as primary (like medaka)
            gene_info['igdb_id'] = primary_id
            gene_info['ens_id'] = secondary_id if secondary_type == 'ens_id' else 'N/A'
        elif primary_type == 'ens_id':
            # Standard species with ENS as primary
            gene_info['ens_id'] = primary_id
            gene_info['igdb_id'] = secondary_id if secondary_type == 'igdb_id' else 'N/A'
        else:
            # Unknown or other format - use generic mapping
            gene_info['ens_id'] = primary_id if primary_type == 'ens_id' else (secondary_id if secondary_type == 'ens_id' else 'N/A')
            gene_info['igdb_id'] = primary_id if primary_type == 'igdb_id' else (secondary_id if secondary_type == 'igdb_id' else 'N/A')
