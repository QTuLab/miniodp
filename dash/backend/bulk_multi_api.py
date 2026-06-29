"""
BulkMulti API - Provides Bulk Multi-omics linkage data access
"""

import logging
import os
import sqlite3
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

import pandas as pd

logger = logging.getLogger(__name__)

class BulkMultiAPI:
    def __init__(self, species_key: str, dataset: Optional[str] = None):
        if not species_key:
            raise ValueError("Species key is required and cannot be empty")
        if not isinstance(species_key, str):
            raise TypeError("Species key must be a string")
        
        self.species = species_key.strip()
        if not self.species:
            raise ValueError("Species key cannot be empty or whitespace")
            
        # Use global DASH_DATA_PATH or fall back to legacy path
        data_root = os.environ.get('DASH_DATA_PATH', Path(__file__).parent.parent / "data")
        self.data_dir = Path(data_root) / self.species
        
        self.logger = logger
        
        # Set dataset path within the BulkMulti directory
        self.dataset = dataset.strip() if isinstance(dataset, str) else None
        self._set_dataset_path(self.dataset)
        self._validate_database()
    
    def _set_dataset_path(self, dataset: Optional[str] = None) -> None:
        """Set the dataset path dynamically"""
        if dataset:
            # Use specified dataset in BulkMulti directory
            self.db_path = self.data_dir / "BulkMulti" / dataset / "BulkMulti.db"
            self.atac_signals_dir = self.data_dir / "BulkMulti" / dataset / "atac_signals"
        else:
            # Auto-discover the newest available dataset
            bulkmulti_dir = self.data_dir / "BulkMulti"
            if bulkmulti_dir.exists():
                available_datasets = [d.name for d in bulkmulti_dir.iterdir() if d.is_dir()]
                if available_datasets:
                    # Use newest dataset (first in reverse sorted order)
                    newest_dataset = sorted(available_datasets, reverse=True)[0]
                    self.dataset = newest_dataset
                    self.db_path = bulkmulti_dir / newest_dataset / "BulkMulti.db"
                    self.atac_signals_dir = bulkmulti_dir / newest_dataset / "atac_signals"
                    self.logger.info(f"Auto-selected dataset: {newest_dataset}")
                else:
                    self.logger.warning(f"No datasets found in {bulkmulti_dir}")
                    self.db_path = None
                    self.atac_signals_dir = None
            else:
                self.logger.warning(f"BulkMulti directory not found: {bulkmulti_dir}")
                self.db_path = None
                self.atac_signals_dir = None
    
    def set_dataset(self, dataset: Optional[str]) -> None:
        """Change the active dataset"""
        if dataset is not None and not isinstance(dataset, str):
            raise TypeError("Dataset name must be a string or None")
        normalized = dataset.strip() if isinstance(dataset, str) else None
        self.dataset = normalized or None
        self._set_dataset_path(self.dataset)
        self._validate_database()
    
    def _assert_data_available(self):
        if not self.db_path or not self.db_path.exists():
            raise FileNotFoundError(f"BulkMulti dataset not found for species {self.species}. Expected database at: {self.db_path}")

    def _connect_ro(self):
        if not self.db_path:
            raise FileNotFoundError("BulkMulti database path is not set")
        return sqlite3.connect(f"file:{self.db_path}?mode=ro", uri=True, check_same_thread=False)
        
    def _validate_database(self):
        """Validate that BulkMulti.db exists and has required tables"""
        if not self.db_path or not self.db_path.exists():
            self.logger.warning(f"BulkMulti database not found: {self.db_path}")
            return
            
        try:
            with self._connect_ro() as conn:
                tables = pd.read_sql_query("SELECT name FROM sqlite_master WHERE type='table'", conn)
                required_tables = ['linkages', 'genes', 'exons', 'peaks', 'rna_expression']
                
                missing_tables = [t for t in required_tables if t not in tables['name'].values]
                if missing_tables:
                    self.logger.warning(f"Missing tables in database: {missing_tables}")
                else:
                    self.logger.info(f"✅ BulkMulti database validated: {len(required_tables)} tables found")
        except Exception as e:
            self.logger.error(f"Error validating database: {e}")
    
    def _normalize_gene_names(self, gene_names: Optional[Sequence[str]]) -> List[str]:
        """Validate and normalize gene name inputs."""
        if gene_names is None:
            return []
        if isinstance(gene_names, str):
            candidates = [gene_names]
        elif not isinstance(gene_names, Sequence):
            raise TypeError("gene_names must be a sequence of strings")
        else:
            candidates = list(gene_names)
        
        normalized: List[str] = []
        for raw in candidates:
            if raw is None:
                continue
            if not isinstance(raw, str):
                raise TypeError(f"Gene name must be a string, got {type(raw)!r}")
            stripped = raw.strip()
            if stripped:
                normalized.append(stripped)
        return normalized

    def _normalize_value_threshold(self, value_threshold: float) -> float:
        """Ensure linkage threshold is a non-negative float."""
        if isinstance(value_threshold, bool):
            raise TypeError("value_threshold cannot be a boolean")
        try:
            threshold = float(value_threshold)
        except (TypeError, ValueError) as exc:
            raise TypeError("value_threshold must be numeric") from exc
        if threshold < 0:
            raise ValueError("value_threshold must be greater than or equal to zero")
        return threshold

    def _validate_region_bounds(self, start: Optional[int], end: Optional[int]) -> Tuple[Optional[int], Optional[int]]:
        """Validate genomic region bounds when optional coordinates are provided."""
        def _coerce(value: Optional[int], label: str) -> Optional[int]:
            if value is None:
                return None
            if isinstance(value, bool):
                raise TypeError(f"{label} must not be boolean")
            try:
                coerced = int(value)
            except (TypeError, ValueError) as exc:
                raise TypeError(f"{label} must be an integer") from exc
            if coerced < 0:
                raise ValueError(f"{label} must be non-negative")
            return coerced
        
        start_val = _coerce(start, "start")
        end_val = _coerce(end, "end")
        if start_val is not None and end_val is not None and start_val > end_val:
            raise ValueError("start cannot be greater than end")
        return start_val, end_val

    def _require_chromosome(self, chromosome: Optional[str]) -> Optional[str]:
        """Ensure chromosome string is valid when required."""
        if chromosome is None:
            return None
        if not isinstance(chromosome, str):
            raise TypeError("chromosome must be a string when provided")
        cleaned = chromosome.strip()
        if not cleaned:
            raise ValueError("chromosome cannot be empty when provided")
        return cleaned
    
    def get_gene_linkages(
        self,
        gene_names: Optional[Sequence[str]],
        value_threshold: float = 0.1
    ) -> pd.DataFrame:
        """Get peak-to-gene linkages for specific genes"""
        try:
            self._assert_data_available()
            normalized_genes = self._normalize_gene_names(gene_names)
            if not normalized_genes:
                self.logger.debug("No gene names provided for linkage query; returning empty DataFrame")
                return pd.DataFrame()
            threshold = self._normalize_value_threshold(value_threshold)
            
            with self._connect_ro() as conn:
                placeholders = ','.join(['?' for _ in normalized_genes])
                query = f"""
                    SELECT seqnames as Chromosome, start as Start, end as End, 
                           gene_name as Gene_name, value as Value, dist as Dist
                    FROM linkages 
                    WHERE gene_name IN ({placeholders}) 
                    AND value >= ?
                    ORDER BY value DESC
                """
                
                gene_linkages = pd.read_sql_query(
                    query, conn, 
                    params=normalized_genes + [threshold]
                )
            
            self.logger.info("Found %d linkages for %d genes", len(gene_linkages), len(normalized_genes))
            return gene_linkages
            
        except Exception as e:
            self.logger.error(f"Error getting gene linkages: {e}")
            raise

    def get_gene_linkages_by_gene(
        self,
        gene_names: Optional[Sequence[str]],
        value_threshold: float = 0.1
    ) -> Dict[str, pd.DataFrame]:
        """Return linkage DataFrames keyed by gene name using a single SQL query."""
        linkages = self.get_gene_linkages(gene_names, value_threshold)
        if linkages.empty:
            return {}

        grouped: Dict[str, pd.DataFrame] = {}
        for gene_name, subset in linkages.groupby('Gene_name'):
            grouped[gene_name] = subset.reset_index(drop=True)

        return grouped
    
    def get_gene_info(self, gene_names: Optional[Sequence[str]]) -> pd.DataFrame:
        """Get gene information for specific genes"""
        try:
            self._assert_data_available()
            normalized_genes = self._normalize_gene_names(gene_names)
            if not normalized_genes:
                self.logger.debug("No gene names provided for gene info query; returning empty DataFrame")
                return pd.DataFrame()
            
            with self._connect_ro() as conn:
                placeholders = ','.join(['?' for _ in normalized_genes])
                query = f"""
                    SELECT * FROM genes 
                    WHERE Symbol IN ({placeholders})
                """
                
                gene_info = pd.read_sql_query(query, conn, params=normalized_genes)
            
            self.logger.info("Found %d genes from %d requested", len(gene_info), len(normalized_genes))
            return gene_info
            
        except Exception as e:
            self.logger.error(f"Error getting gene info: {e}")
            raise
    
    def get_exon_data(
        self,
        gene_names: Optional[Sequence[str]],
        chromosome: Optional[str] = None,
        start: Optional[int] = None,
        end: Optional[int] = None
    ) -> pd.DataFrame:
        """Get real exon coordinates for genes"""
        try:
            self._assert_data_available()
            normalized_genes = self._normalize_gene_names(gene_names)
            chromosome_value = self._require_chromosome(chromosome)
            start_val, end_val = self._validate_region_bounds(start, end)
            
            conditions = []
            params = []
            
            if normalized_genes:
                placeholders = ','.join(['?' for _ in normalized_genes])
                conditions.append(f"Gene_Symbol IN ({placeholders})")
                params.extend(normalized_genes)
                
            if chromosome_value:
                conditions.append("Chromosome = ?")
                params.append(chromosome_value)
                
            if start_val is not None:
                conditions.append("End >= ?")
                params.append(start_val)
                
            if end_val is not None:
                conditions.append("Start <= ?")
                params.append(end_val)
                
            where_clause = " AND ".join(conditions) if conditions else "1=1"
            
            query = f"""
                SELECT * FROM exons 
                WHERE {where_clause}
                ORDER BY Chromosome, Start
            """
            
            self.logger.debug("🔍 Exon Query params: %s", params)
            self.logger.debug("🔍 Exon Query param types: %s", [type(p) for p in params])
            clean_params = []
            for p in params:
                if isinstance(p, str):
                    clean_params.append(str(p))
                else:
                    clean_params.append(int(p))
            
            with self._connect_ro() as conn:
                exon_data = pd.read_sql_query(query, conn, params=clean_params)
            
            self.logger.info(f"Found {len(exon_data)} exons for region query")
            if not exon_data.empty:
                self.logger.debug("🔍 Exon columns: %s", exon_data.columns.tolist())
                sample = exon_data.iloc[0].to_dict() if len(exon_data) > 0 else 'None'
                self.logger.debug("🔍 Sample exon row: %s", sample)
            return exon_data
            
        except Exception as e:
            self.logger.error(f"Error getting exon data: {e}")
            return pd.DataFrame()
    
    def get_gene_structure(self, gene_name: str) -> Dict:
        """Get complete gene structure including all exons"""
        try:
            normalized = self._normalize_gene_names([gene_name])
            if not normalized:
                raise ValueError("gene_name is required")
            gene_info = self.get_gene_info(normalized)
            exon_data = self.get_exon_data(normalized)
            
            if gene_info.empty:
                return {}
                
            gene = gene_info.iloc[0]
            structure = {
                'gene': {
                    'symbol': gene['Symbol'],
                    'chromosome': gene['Chromosome'],
                    'start': gene['Start'],
                    'end': gene['End'],
                    'strand': gene['Strand']
                },
                'exons': []
            }
            
            for _, exon in exon_data.iterrows():
                structure['exons'].append({
                    'start': exon['Start'],
                    'end': exon['End'],
                    'strand': exon['Strand']
                })
                
            return structure
            
        except Exception as e:
            self.logger.error(f"Error getting gene structure: {e}")
            return {}
    
    def get_genomic_region_data(
        self,
        chromosome: str,
        start: int,
        end: int,
        gene_names: Optional[Sequence[str]]
    ) -> Dict:
        """Get genomic region data for visualization"""
        try:
            self._assert_data_available()
            chromosome_value = self._require_chromosome(chromosome)
            start_val, end_val = self._validate_region_bounds(start, end)
            if start_val is None or end_val is None:
                raise ValueError("start and end positions are required")
            normalized_genes = self._normalize_gene_names(gene_names)
            
            linkages = self.get_gene_linkages(normalized_genes)
            gene_info = self.get_gene_info(normalized_genes)
            
            region_linkages = linkages[
                (linkages['Chromosome'] == chromosome_value) &
                (linkages['Start'] <= end_val) &
                (linkages['End'] >= start_val)
            ]
            region_genes = gene_info[
                (gene_info['Chromosome'] == chromosome_value) &
                (gene_info['Start'] <= end_val) &
                (gene_info['End'] >= start_val)
            ]
            
            region_features = pd.DataFrame()
            try:
                with self._connect_ro() as conn:
                    region_features = pd.read_sql_query(
                        """
                        SELECT seqnames as Chromosome, start as Start, end as End
                        FROM peaks
                        WHERE seqnames = ?
                        AND start <= ?
                        AND end >= ?
                        """,
                        conn,
                        params=[chromosome_value, end_val, start_val]
                    )
            except Exception as feature_err:
                self.logger.warning(f"Error querying peaks table for region data: {feature_err}")
            
            return {
                'linkages': region_linkages,
                'genes': region_genes,
                'features': region_features,
                'chromosome': chromosome_value,
                'start': start_val,
                'end': end_val
            }
            
        except Exception as e:
            self.logger.error(f"Error getting genomic region data: {e}")
            return {}
    
    def get_linkage_statistics(self, gene_names: Optional[Sequence[str]]) -> Dict:
        """Get statistics for peak-to-gene linkages"""
        try:
            normalized_genes = self._normalize_gene_names(gene_names)
            if not normalized_genes:
                return {}
            linkages = self.get_gene_linkages(normalized_genes, value_threshold=0.0)
            
            if linkages.empty:
                return {}
            
            stats = {}
            for gene in normalized_genes:
                gene_data = linkages[linkages['Gene_name'] == gene]
                if not gene_data.empty:
                    stats[gene] = {
                        'total_linkages': len(gene_data),
                        'strong_linkages': len(gene_data[gene_data['Value'] >= 0.5]),
                        'moderate_linkages': len(gene_data[(gene_data['Value'] >= 0.2) & (gene_data['Value'] < 0.5)]),
                        'weak_linkages': len(gene_data[gene_data['Value'] < 0.2]),
                        'max_linkage_value': float(gene_data['Value'].max()),
                        'mean_linkage_value': float(gene_data['Value'].mean()),
                        'median_distance': float(gene_data['Dist'].median()),
                        'chromosomes': gene_data['Chromosome'].unique().tolist()
                    }
            
            return stats
            
        except Exception as e:
            self.logger.error(f"Error getting linkage statistics: {e}")
            return {}
    
    def get_chromosome_view_data(
        self,
        gene_names: Optional[Sequence[str]],
        chromosome: Optional[str] = None
    ) -> Dict:
        """Get chromosome-level view data for genes"""
        try:
            normalized_genes = self._normalize_gene_names(gene_names)
            if not normalized_genes:
                return {}
            linkages = self.get_gene_linkages(normalized_genes)
            gene_info = self.get_gene_info(normalized_genes)
            
            if linkages.empty or gene_info.empty:
                return {}
            
            # If chromosome not specified, use the first gene's chromosome
            if chromosome is None and not gene_info.empty:
                chromosome = gene_info.iloc[0]['Chromosome']
            chromosome_value = self._require_chromosome(chromosome)
            
            # Filter by chromosome
            chr_linkages = linkages[linkages['Chromosome'] == chromosome_value]
            chr_genes = gene_info[gene_info['Chromosome'] == chromosome_value]
            
            # Calculate region bounds
            if not chr_genes.empty:
                gene_starts = chr_genes['Start'].min()
                gene_ends = chr_genes['End'].max()
                
                # Expand region to include linked peaks
                if not chr_linkages.empty:
                    peak_starts = chr_linkages['Start'].min()
                    peak_ends = chr_linkages['End'].max()
                    region_start = min(gene_starts, peak_starts) - 50000
                    region_end = max(gene_ends, peak_ends) + 50000
                else:
                    region_start = gene_starts - 50000
                    region_end = gene_ends + 50000
                
                return {
                    'chromosome': chromosome_value,
                    'region_start': int(max(0, region_start)),
                    'region_end': int(region_end),
                    'genes': chr_genes,
                    'linkages': chr_linkages,
                    'gene_count': len(chr_genes),
                    'linkage_count': len(chr_linkages)
                }
            
            return {}
            
        except Exception as e:
            self.logger.error(f"Error getting chromosome view data: {e}")
            return {}
    
    def search_genes_by_region(self, chromosome: str, start: int, end: int) -> List[str]:
        """Search for genes in a genomic region"""
        try:
            self._assert_data_available()
            chromosome_value = self._require_chromosome(chromosome)
            start_val, end_val = self._validate_region_bounds(start, end)
            if start_val is None or end_val is None:
                raise ValueError("start and end positions are required")
            
            query = """
                SELECT Symbol FROM genes 
                WHERE Chromosome = ? 
                AND Start <= ? 
                AND End >= ?
            """
            
            self.logger.info(f"🔍 Searching genes in region: {chromosome_value}:{start_val}-{end_val}")
            
            params = [chromosome_value, end_val, start_val]
            self.logger.info(f"🔍 Query params: {params}")
            self.logger.info(f"🔍 Query param types: {[type(p) for p in params]}")
            
            with self._connect_ro() as conn:
                region_genes = pd.read_sql_query(
                    query, conn, 
                    params=params
                )
            
            gene_list = region_genes['Symbol'].tolist()
            self.logger.info(f"🔍 Found {len(gene_list)} genes in region: {gene_list}")
            
            return gene_list
            
        except Exception as e:
            self.logger.error(f"Error searching genes by region: {e}")
            return []
    
    def get_available_chromosomes(self) -> List[str]:
        """Get list of available chromosomes"""
        try:
            self._assert_data_available()
            
            with self._connect_ro() as conn:
                chromosomes = pd.read_sql_query(
                    "SELECT DISTINCT seqnames as Chromosome FROM linkages ORDER BY seqnames", 
                    conn
                )
            
            return chromosomes['Chromosome'].tolist()
            
        except Exception as e:
            self.logger.error(f"Error getting available chromosomes: {e}")
            return []

    def get_global_statistics(self) -> Dict:
        """Get global Peak2Gene statistics"""
        try:
            self._assert_data_available()
            
            stats = {}
            with self._connect_ro() as conn:
                linkages_df = pd.read_sql_query(
                    "SELECT gene_name as Gene_name, seqnames as Chromosome, value as Value, dist as Dist FROM linkages",
                    conn
                )
                genes_df = pd.read_sql_query(
                    "SELECT Symbol, Chromosome FROM genes",
                    conn
                )
                peaks_df = pd.read_sql_query(
                    "SELECT seqnames as Chromosome FROM peaks",
                    conn
                )
            
            if not linkages_df.empty:
                stats['linkages'] = {
                    'total_linkages': len(linkages_df),
                    'unique_genes': linkages_df['Gene_name'].nunique(),
                    'chromosomes': linkages_df['Chromosome'].nunique(),
                    'value_range': [
                        float(linkages_df['Value'].min()),
                        float(linkages_df['Value'].max())
                    ],
                    'distance_range': [
                        int(linkages_df['Dist'].min()),
                        int(linkages_df['Dist'].max())
                    ],
                    'strong_linkages': int((linkages_df['Value'] >= 0.5).sum()),
                    'moderate_linkages': int(((linkages_df['Value'] >= 0.2) & (linkages_df['Value'] < 0.5)).sum()),
                    'weak_linkages': int((linkages_df['Value'] < 0.2).sum())
                }
            
            if not genes_df.empty:
                stats['genes'] = {
                    'total_genes': len(genes_df),
                    'chromosomes': genes_df['Chromosome'].nunique()
                }
            
            if not peaks_df.empty:
                stats['features'] = {
                    'total_features': len(peaks_df),
                    'chromosomes': peaks_df['Chromosome'].nunique()
                }
            
            return stats
            
        except Exception as e:
            self.logger.error(f"Error getting global statistics: {e}")
            return {}
