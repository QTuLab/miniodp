"""
Expression Data API - Parquet-based implementation for better performance
"""

import pandas as pd
import numpy as np
from typing import List, Dict, Optional, Tuple
from pathlib import Path
import logging
import time


logger = logging.getLogger(__name__)

class ExpressionAPI:
    def __init__(self, species_key: str, adapter=None):
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
        self.parquet_dir = self.data_dir / "bulkRNA"
        self.data_file = self.parquet_dir / "expression_data.parquet"
        self.gene_index_file = self.parquet_dir / "gene_index.parquet"
        self.source_index_file = self.parquet_dir / "source_index.parquet"
        
        if not self.data_file.exists():
            raise FileNotFoundError(f"Parquet data not found: {self.data_file}")
        
        # Set up species adapter for ID handling
        if adapter:
            self.adapter = adapter
        else:
            from .species_adapters import create_species_adapter
            self.adapter = create_species_adapter(species_key)
        
        # Get the gene ID column name used in data storage (common_id column)
        self.gene_id_column = self.adapter.get_common_id_column()
        
        self.logger = logger
        
        # Load indexes into memory (small files)
        self.logger.info("📊 Loading indexes...")
        self.gene_index = pd.read_parquet(self.gene_index_file)
        self.source_index = pd.read_parquet(self.source_index_file)
        self.logger.info("✅ Loaded %d genes, %d sources", len(self.gene_index), len(self.source_index))
        
        try:
            samples_df = pd.read_parquet(self.data_file, columns=["Source", "Sample"])
            samples_df["Sample_Source"] = samples_df["Source"] + "::" + samples_df["Sample"]
            self.samples_index = samples_df.drop_duplicates()
        except Exception as exc:
            self.logger.warning("Failed to preload samples index: %s", exc)
            self.samples_index = pd.DataFrame(columns=["Source", "Sample", "Sample_Source"])

    def get_available_sources(self) -> List[Dict[str, any]]:
        """Get all available studies with metadata"""
        try:
            sources = []
            for source_name, row in self.source_index.iterrows():
                sources.append({
                    'source': source_name,
                    'sample_count': int(row['samples']),
                    'gene_count': int(row['genes']),
                    'total_records': int(row['records']),
                    'mean_tpm': float(row['mean_tpm'])
                })
            
            # Sort by name (descending, so year-prefixed names show newest first)
            sources.sort(key=lambda x: x['source'], reverse=True)
            return sources
            
        except Exception as e:
            self.logger.error(f"Error getting sources: {e}")
            return []

    def get_samples_by_sources(self, sources: List[str]) -> List[Dict[str, any]]:
        """Get samples for specified sources"""
        try:
            df = self.samples_index
            if sources:
                df = df[df["Source"].isin(sources)]
            samples = [
                {
                    "source": row["Source"],
                    "sample": row["Sample"],
                    "note": f"{row['Sample']} ({row['Source']})",
                }
                for _, row in df.iterrows()
            ]
            return samples
            
        except Exception as e:
            self.logger.error(f"Error getting samples: {e}")
            return []

    def get_expression_matrix(self, 
                            gene_ids: List[str], 
                            sources: Optional[List[str]] = None,
                            samples: Optional[List[str]] = None,
                            min_tpm: float = 0.0) -> Tuple[pd.DataFrame, Dict[str, any]]:
        """Get expression matrix for specified genes"""
        if not gene_ids:
            return pd.DataFrame(), {'error': 'No genes provided'}

        try:
            self.logger.info("Loading data for %d genes", len(gene_ids))
            start_time = time.time()
            
            # Build filters for efficient parquet reading using species-specific gene ID column
            filters = [(self.gene_id_column, 'in', gene_ids)]
            if sources:
                filters.append(('Source', 'in', sources))
            
            # Read filtered data from parquet
            df = pd.read_parquet(self.data_file, filters=filters)
            
            use_sample_source_filter = bool(samples and any('::' in sample for sample in samples))
            
            if samples:
                if use_sample_source_filter:
                    if 'Sample_Source' not in df.columns:
                        self.logger.warning("Sample_Source column missing, creating on-the-fly for filtering")
                        df['Sample_Source'] = df['Source'] + '::' + df['Sample']
                    df = df[df['Sample_Source'].isin(samples)]
                else:
                    df = df[df['Sample'].isin(samples)]
            
            if min_tpm > 0:
                df = df[df['TPM'] >= min_tpm]
            
            load_time = time.time() - start_time
            self.logger.info("Loaded %d records in %.3fs", len(df), load_time)
            
            if df.empty:
                return pd.DataFrame(), {'error': 'No data found for query'}
            
            # Use pre-built Sample_Source column (should already exist in data)
            if 'Sample_Source' not in df.columns:
                # Fallback: create if missing (for backward compatibility)
                self.logger.warning("Sample_Source column missing, creating on-the-fly")
                df['Sample_Source'] = df['Source'] + '::' + df['Sample']
            
            # Create expression matrix using the species-specific gene ID column
            matrix = df.pivot_table(
                index=self.gene_id_column,
                columns='Sample_Source',
                values='TPM',
                fill_value=0.0,
                aggfunc='mean'  # Handle duplicates if any
            )
            
            # Calculate metadata
            metadata = {
                'genes_found': len(matrix.index),
                'genes_requested': len(gene_ids),
                'samples_count': len(matrix.columns),
                'total_values': matrix.size,
                'non_zero_values': (matrix > 0).sum().sum(),
                'missing_genes': [gene for gene in gene_ids if gene not in matrix.index],
                'tpm_range': [float(matrix.min().min()), float(matrix.max().max())],
                'sources_included': sources or [],
                'samples_included': samples or [],
                'load_time_seconds': load_time
            }
            
            result = (matrix, metadata)
            return result
            
        except Exception as e:
            self.logger.error(f"Error getting expression matrix: {e}")
            return pd.DataFrame(), {'error': str(e)}

    def transform_expression_data(self, 
                                matrix: pd.DataFrame, 
                                method: str = 'log2',
                                pseudocount: float = 1.0) -> pd.DataFrame:
        """Transform expression data using various methods"""
        if matrix.empty:
            return matrix
        
        if method == 'raw':
            return matrix
        elif method == 'log2':
            return np.log2(matrix + pseudocount)
        elif method == 'zscore_rows':
            return matrix.subtract(matrix.mean(axis=1), axis=0).div(matrix.std(axis=1), axis=0).fillna(0)
        elif method == 'zscore_cols':
            return matrix.subtract(matrix.mean(axis=0), axis=1).div(matrix.std(axis=0), axis=1).fillna(0)
        else:
            return matrix
