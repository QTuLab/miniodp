"""
Abstract base class for species-specific adapters.

This defines the interface that all species adapters must implement,
ensuring consistent behavior across different species while allowing
for species-specific customizations.
"""

import logging
import os
import re
import sqlite3
from abc import ABC, abstractmethod
from pathlib import Path
from typing import Any, Dict, List, Optional

import pandas as pd

try:
    import tomllib
except ModuleNotFoundError:  # Python < 3.11
    import tomli as tomllib


logger = logging.getLogger(__name__)


class SpeciesAdapter(ABC):
    """Abstract base class for species-specific data adapters."""
    
    def __init__(self, species_key: str):
        self.species_key = species_key
        self.logger = logging.getLogger(f"{self.__class__.__module__}.{self.__class__.__name__}.{species_key}")
        self.config = self._load_configuration()
        self._search_field_configs = self._build_search_field_configs()
    
    @abstractmethod
    def get_search_fields(self) -> List[Dict[str, str]]:
        """
        Return available search field configurations for this species.
        
        Returns:
            List of dicts with 'label' and 'value' keys for dropdown options
        """
        pass
    
    @abstractmethod
    def search_genes(self, db_conn: sqlite3.Connection, query: str, search_type: str) -> pd.DataFrame:
        """
        Execute gene search and return standardized results.
        
        Args:
            db_conn: Database connection object
            query: Search query string
            search_type: Type of search to perform
            
        Returns:
            DataFrame with standardized gene information
        """
        pass
    
    @abstractmethod
    def batch_search_genes(self, db_conn: sqlite3.Connection, gene_ids: List[str], search_type: str) -> pd.DataFrame:
        """
        Execute gene search for multiple genes and return standardized results.
        
        Args:
            db_conn: Database connection object
            gene_ids: List of gene identifiers to search for
            search_type: Type of search to perform (e.g., 'ens_id', 'igdb_id')
            
        Returns:
            DataFrame with standardized gene information
        """
        pass
    
    @abstractmethod
    def get_primary_key_column(self) -> str:
        """
        Return the column name used as primary key for expression data.
        
        Returns:
            Column name (e.g., 'ENS' for standard species, 'IGDB' for medaka)
        """
        pass
    
    @abstractmethod
    def get_display_columns(self) -> List[str]:
        """
        Return list of columns to display in results table.
        
        Returns:
            List of column names in display order
        """
        pass
    
    @abstractmethod
    def format_gene_label(self, row: pd.Series) -> str:
        """
        Generate display label for a gene.
        
        Args:
            row: DataFrame row containing gene information
            
        Returns:
            Formatted gene label string
        """
        pass
    
    @abstractmethod
    def get_common_id_column(self) -> str:
        """
        Return the column name for common_id in data storage (parquet files).

        common_id is the species' commonly used ID system.
        For standard species this is 'ENS', for medaka this is 'IGDB'.

        Returns:
            Column name used in parquet files (e.g., 'ENS', 'IGDB')
        """
        pass

    @abstractmethod
    def parse_primary_id(self, primary_id: str) -> Dict[str, str]:
        """
        Parse primary_id into component IDs, extract common_id.

        Args:
            primary_id: The primary ID (e.g., 'IGDB::ENS' for medaka, 'ENS' for standard)

        Returns:
            Dict with 'igdb' and 'ens' keys (some may be empty strings)
        """
        pass

    @abstractmethod
    def get_common_ids(self, primary_ids: List[str]) -> List[str]:
        """
        Extract common_ids from primary_ids for data queries.

        This converts primary_ids to the IDs used in parquet files.
        For medaka: extracts IGDB from 'IGDB::ENS' format
        For standard: returns ENS IDs unchanged

        Args:
            primary_ids: List of primary IDs from the shopping cart

        Returns:
            List of common_ids suitable for expression data queries
        """
        pass

    @abstractmethod
    def get_gene_sequences(self, db_conn: sqlite3.Connection, gene_id: str) -> List[Dict[str, str]]:
        """
        Get all transcript sequences for a gene.
        
        Args:
            db_conn: Database connection object
            gene_id: Gene identifier to get sequences for
            
        Returns:
            List of dicts with 'transcript_id', 'length', and 'fasta' keys
        """
        pass
    
    def validate_search_type(self, search_type: str) -> bool:
        """
        Validate if search type is supported by this adapter.
        
        Args:
            search_type: Type of search to validate
            
        Returns:
            True if supported, False otherwise
        """
        valid_types = [field['value'] for field in self.get_search_fields()]
        return search_type in valid_types
    
    def get_species_info(self) -> Dict[str, Any]:
        """
        Return basic information about this species.
        
        Returns:
            Dict with species metadata
        """
        return {
            'species_key': self.species_key,
            'adapter_type': self.__class__.__name__,
            'primary_key': self.get_primary_key_column(),
            'search_fields': len(self.get_search_fields()),
            'display_columns': len(self.get_display_columns())
        }
    
    @abstractmethod
    def format_search_results(self, raw_results: pd.DataFrame, db_conn: sqlite3.Connection = None) -> List[Dict[str, Any]]:
        """
        Format raw search results into standardized format.

        Args:
            raw_results: Raw DataFrame from search_genes()
            db_conn: Optional database connection for additional lookups (e.g., IGDB_NAME)

        Returns:
            List of standardized result dictionaries with keys:
            - primary_id: Main identifier for this species
            - secondary_id: Optional secondary identifier
            - gene_name: Gene name
            - description: Gene description
            - gene_type: Gene type
            - locus: Genomic location (optional)
        """
        pass
    
    @abstractmethod
    def get_external_links(self, result: Dict[str, Any]) -> Dict[str, str]:
        """
        Generate external database links for a search result.
        
        Args:
            result: Standardized result dictionary
            
        Returns:
            Dict mapping link names to URLs
        """
        pass
    
    @abstractmethod
    def get_tf_families(self, db_conn: sqlite3.Connection, primary_id: str, secondary_id: Optional[str] = None) -> List[str]:
        """
        Get transcription factor families for a gene.
        
        Args:
            db_conn: Database connection
            primary_id: Primary gene identifier
            secondary_id: Secondary gene identifier (optional)
            
        Returns:
            List of TF family names
        """
        pass
    
    @abstractmethod
    def get_summary_columns(self) -> Dict[str, str]:
        """
        Return summary table column configuration for this species.
        
        Returns:
            Dict mapping column keys to display names
        """
        pass
    
    @abstractmethod
    def format_search_result_display(self, result: Dict[str, Any]) -> Dict[str, Any]:
        """
        Format search result for display in search results list.
        
        Args:
            result: Search result dictionary with gene information
            
        Returns:
            Dict with display configuration:
            - main_text: Primary display text (gene name)
            - secondary_items: List of secondary info items to display
        """
        pass
    
    @abstractmethod
    def is_tf_search_field(self, field: str) -> bool:
        """
        Check if a search field is for transcription factor searches.
        
        Args:
            field: Search field name to check
            
        Returns:
            True if the field is for TF searches, False otherwise
        """
        pass
    
    @abstractmethod
    def convert_ids_to_gene_names(self, gene_ids: List[str]) -> List[str]:
        """
        Convert gene IDs to gene names for this species.
        
        Args:
            gene_ids: List of gene identifiers (can be any format)
            
        Returns:
            List of gene names corresponding to the input IDs
        """
        pass
    
    @abstractmethod
    def get_column_display_config(self) -> Dict[str, str]:
        """
        Return column display configuration mapping internal field names to display labels.
        
        Returns:
            Dict mapping field names to display labels (e.g., 'primary_id': 'Gene ID')
        """
        pass
    
    def _load_configuration(self) -> Dict[str, Any]:
        """
        Load species configuration from TOML file.
        
        Returns:
            Configuration dictionary, empty if file not found
        """
        try:
            candidates = []

            data_root = os.environ.get("DASH_DATA_PATH")
            if data_root:
                candidates.append(Path(data_root) / self.species_key / "geneinfo.toml")

            base_dir = Path(__file__).resolve().parents[2]
            candidates.append(base_dir / "data" / self.species_key / "geneinfo.toml")

            for toml_path in candidates:
                if toml_path.exists():
                    with open(toml_path, "rb") as f:
                        return tomllib.load(f)

            self.logger.warning(
                "No TOML config found for %s in %s, using adapter defaults",
                self.species_key,
                candidates,
            )
            return {}
        except Exception as e:
            self.logger.error("Error loading config for %s: %s", self.species_key, e)
            return {}
    
    def _load_species_global_config(self) -> Dict[str, Any]:
        """
        Load global species configuration from species_config.toml.
        Merges default settings with species-specific overrides.
        
        Returns:
            Merged species configuration dictionary
        """
        try:
            config_path = Path(__file__).parent.parent.parent / "data" / "species_adapters.toml"
            
            if config_path.exists():
                with open(config_path, "rb") as f:
                    config = tomllib.load(f)
                    
                    # Start with default configuration
                    merged_config = config.get('default', {}).copy()
                    
                    # Override with species-specific settings
                    species_config = config.get('species', {}).get(self.species_key, {})
                    merged_config.update(species_config)
                    
                    return merged_config
            self.logger.warning("Global species config not found: %s", config_path)
            return {}
        except Exception as e:
            self.logger.error("Error loading global species config for %s: %s", self.species_key, e)
            return {}
    
    def _load_hugo_species_config(self) -> Dict[str, Any]:
        """
        Load species configuration from Hugo's data/species_display.toml file.
        
        Returns:
            Species configuration dictionary, empty if not found
        """
        try:
            candidates = []

            # 1) Data root (synced into dash/data/, mounted to container)
            data_root = os.environ.get("DASH_DATA_PATH")
            if data_root:
                candidates.append(Path(data_root) / "species_display.toml")

            # 2) Repo root layout: <repo>/hugo/data/species_display.toml
            repo_root = Path(__file__).resolve().parents[3]
            candidates.append(repo_root / "hugo" / "data" / "species_display.toml")

            # 3) Fallback: sibling "hugo/data" relative to current working dir (runtime override)
            candidates.append(Path.cwd() / "hugo" / "data" / "species_display.toml")

            for path in candidates:
                if path.exists():
                    with open(path, "rb") as f:
                        all_species = tomllib.load(f)
                        return all_species.get(self.species_key, {})

            self.logger.warning("Hugo species config not found in candidates: %s", candidates)
            return {}
        except Exception as e:
            self.logger.error("Error loading Hugo config for %s: %s", self.species_key, e)
            return {}
    
    def _get_ensembl_species_name(self) -> str:
        """
        Get Ensembl species name from Hugo configuration.
        
        Returns:
            Ensembl species name (scientific name with spaces replaced by underscores)
        """
        hugo_config = self._load_hugo_species_config()
        scientific_name = hugo_config.get('scientific_name', '')
        
        if scientific_name:
            # Convert "Danio rerio" -> "Danio_rerio" for Ensembl URLs
            return scientific_name.replace(' ', '_')
        else:
            # Fallback to species key if scientific name not found
            self.logger.warning("No scientific name found for %s, using species key", self.species_key)
            return self.species_key
    
    def _build_search_field_configs(self) -> Dict[str, Dict[str, str]]:
        """
        Build mapping of search field values to their configurations.
        
        Returns:
            Dict mapping field values to config dicts
        """
        if not self.config.get('search_fields'):
            return {}
        
        return {field['value']: field for field in self.config['search_fields']}
    
    def get_available_search_fields(self) -> List[Dict[str, str]]:
        """
        Get available search fields, prioritizing TOML config over adapter defaults.
        
        Returns:
            List of search field configurations
        """
        if self.config.get('search_fields'):
            return [{'label': f['label'], 'value': f['value']} 
                   for f in self.config['search_fields']]
        else:
            return self.get_search_fields()
    
    def _clean_query(self, query: str) -> str:
        """Clean and normalize search query."""
        query = re.sub(r'\s+', ' ', query.strip())
        return re.sub(r'[^\w\s.:-]', '', query)
    
    def _format_locus_string(self, row: pd.Series) -> Optional[str]:
        """Format genomic locus information from row data."""
        if 'chromosome' not in row or pd.isna(row.get('chromosome')):
            return None
            
        chrom = str(row['chromosome'])
        start = row.get('start')
        end = row.get('end')
        strand = row.get('strand')

        if chrom.lower().startswith("chr"):
            display_chrom = chrom
        elif re.fullmatch(r"(?:[0-9]+|[0-9]+[A-Za-z]?|X|Y|MT|M)", chrom):
            display_chrom = f"chr{chrom}"
        else:
            display_chrom = chrom
        
        if start and end:
            locus = f"{display_chrom}:{start}-{end}"
            if strand in [1, '1', '+']:
                locus += " (+)"
            elif strand in [-1, '-1', '-']:
                locus += " (-)"
            return locus
        else:
            return display_chrom
    
    def detect_id_format(self, gene_id: str) -> str:
        """
        Detect the format type of a gene ID using configuration-driven patterns.
        
        Args:
            gene_id: Gene identifier to analyze
            
        Returns:
            Search field type or 'gene_name' if unrecognized
        """
        import re
        
        # Get species configuration
        species_config = self._load_species_global_config()
        id_formats = species_config.get('id_formats', [])
        
        # Check for common patterns first (universal)
        if gene_id.startswith('GO:'):
            return 'go_id'
        elif gene_id.startswith('IPR'):
            return 'interpro_id'
        
        # Check species-specific ID formats
        for id_pattern in id_formats:
            try:
                if re.match(id_pattern, gene_id):
                    # Map pattern to field type based on species configuration
                    primary_key = species_config.get('primary_key', 'ENS')
                    
                    # Determine field type based on pattern matching
                    if 'IGDB' in id_pattern:
                        return 'igdb_id' if primary_key == 'IGDB' else 'secondary_id'
                    elif 'ENS' in id_pattern:
                        return 'ens_id' if primary_key == 'ENS' else 'secondary_id'
                    elif primary_key == 'IGDB':
                        return 'igdb_id'
                    elif primary_key == 'ENS':
                        return 'ens_id'
                    return 'primary_id'
            except re.error:
                # Invalid regex pattern, skip
                continue
        
        # Default to gene name search for unrecognized patterns
        return 'gene_name'
