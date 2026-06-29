"""
Landscape View Visualization for BulkMulti Analysis
Based on R implementation in peak2gene_function.R and plot_function.R
"""

import csv
import logging
import os
import sqlite3
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
import pyBigWig
from plotly.subplots import make_subplots

logger = logging.getLogger(__name__)

class LandscapeViewAPI:
    """
    Landscape View Visualization API for BulkMulti Analysis
    
    This class provides a unified interface for creating multi-track landscape view plots
    with ATAC-seq accessibility data, gene annotations, regulatory elements, and linkage arcs.
    """
    
    # Class-level configuration constants
    DEFAULT_FONT_FAMILY = "Arial"
    DEFAULT_PANEL_TITLE_SIZE = 18
    DEFAULT_TRACK_LABEL_SIZE = 18
    DEFAULT_GENE_LABEL_SIZE  = 16
    
    # Default track heights (in pixels)
    DEFAULT_ATAC_HEIGHT    =  60
    DEFAULT_GENE_HEIGHT    = 200
    DEFAULT_ELS_HEIGHT     =  30
    DEFAULT_LINKAGE_HEIGHT = 300
    DEFAULT_PANEL_TITLE_HEIGHT = 80
    
    # ATAC visualization defaults
    DEFAULT_ATAC_HEIGHT_RATIO = 0.2
    
    def __init__(self, species_key: str, custom_config: dict = None, dataset: str = None):
        if not species_key:
            raise ValueError("Species key is required and cannot be empty")
        if not isinstance(species_key, str):
            raise TypeError("Species key must be a string")
        
        self.species = species_key.strip()
        if not self.species:
            raise ValueError("Species key cannot be empty or whitespace")
            
        self.dataset = dataset
        self.logger = logger
        
        # Use global DASH_DATA_PATH or fall back to legacy path
        data_root = os.environ.get('DASH_DATA_PATH', Path(__file__).parent.parent / "data")
        species_data_dir = Path(data_root) / species_key
        
        # Set data directory - BulkMulti format only
        if dataset:
            self.data_dir = species_data_dir / "BulkMulti" / dataset
            self.atac_signals_dir = self.data_dir / "atac_signals"
        else:
            # Auto-discover the newest available dataset
            bulkmulti_dir = species_data_dir / "BulkMulti"
            if bulkmulti_dir.exists():
                available_datasets = [d.name for d in bulkmulti_dir.iterdir() if d.is_dir()]
                if available_datasets:
                    # Use newest dataset (first in reverse sorted order)
                    newest_dataset = sorted(available_datasets, reverse=True)[0]
                    self.dataset = newest_dataset
                    self.data_dir = bulkmulti_dir / newest_dataset
                    self.atac_signals_dir = self.data_dir / "atac_signals"
                    self.logger.info(f"Auto-selected dataset: {newest_dataset}")
                else:
                    raise ValueError(f"No datasets found in BulkMulti directory: {bulkmulti_dir}")
            else:
                raise FileNotFoundError(f"BulkMulti directory not found: {bulkmulti_dir}")
        
        # Apply custom configuration or use defaults
        config = custom_config or {}
        
        # Font settings
        self.main_font_family = config.get('font_family', self.DEFAULT_FONT_FAMILY)
        self.panel_title_font_size = config.get('panel_title_size', self.DEFAULT_PANEL_TITLE_SIZE)
        self.track_label_font_size = config.get('track_label_size', self.DEFAULT_TRACK_LABEL_SIZE)
        self.gene_label_font_size = config.get('gene_label_size', self.DEFAULT_GENE_LABEL_SIZE)
        
        # ATAC visualization settings
        self.atac_max_height_ratio = config.get('atac_height_ratio', self.DEFAULT_ATAC_HEIGHT_RATIO)
        self.show_els_guide_lines = config.get('show_els_guide_lines', False)
        
        # Track heights (in pixels)
        self.atac_track_height = config.get('atac_height', self.DEFAULT_ATAC_HEIGHT)
        self.gene_track_height = config.get('gene_height', self.DEFAULT_GENE_HEIGHT)
        self.els_track_height = config.get('els_height', self.DEFAULT_ELS_HEIGHT)
        self.linkage_track_height = config.get('linkage_height', self.DEFAULT_LINKAGE_HEIGHT)
        self.panel_title_height = config.get('panel_title_height', self.DEFAULT_PANEL_TITLE_HEIGHT)
        
        # Database paths
        if self.data_dir:
            self.bulkmulti_db_path = self.data_dir / "BulkMulti.db"
        else:
            self.bulkmulti_db_path = None
        
        # Sample order will be loaded from configuration file or auto-detected
        self.sample_order = []
        
        # Initialize sample configuration from file or auto-detection
        self.sample_config = {}
        self._load_sample_configuration()
        
        # Dynamic color generation for variable track counts
        self._generate_sample_colors()

    def _connect_ro(self, db_path: Path):
        return sqlite3.connect(f"file:{db_path}?mode=ro", uri=True, check_same_thread=False)

    def _get_render_settings(self, overrides: Optional[Dict[str, Any]] = None) -> Dict[str, Any]:
        """
        Build a render configuration dictionary without mutating shared state.
        Allows per-request overrides (e.g., slider adjustments) while keeping
        the API instance reusable across callbacks.
        """
        settings: Dict[str, Any] = {
            'font_family': self.main_font_family,
            'panel_title_font_size': self.panel_title_font_size,
            'track_label_font_size': self.track_label_font_size,
            'gene_label_font_size': self.gene_label_font_size,
            'atac_height_ratio': self.atac_max_height_ratio,
            'show_els_guide_lines': self.show_els_guide_lines,
            'atac_track_height': self.atac_track_height,
            'gene_track_height': self.gene_track_height,
            'els_track_height': self.els_track_height,
            'linkage_track_height': self.linkage_track_height,
            'panel_title_height': self.panel_title_height,
        }

        if overrides:
            for key, value in overrides.items():
                if key in settings and value is not None:
                    settings[key] = value

        return settings
    
    def set_dataset(self, dataset: Optional[str]) -> None:
        """Change the active dataset"""
        if dataset is not None and not isinstance(dataset, str):
            raise TypeError("dataset must be a string or None")
        self.dataset = dataset.strip() if isinstance(dataset, str) else dataset
        if not self.dataset:
            raise ValueError("dataset is required when switching BulkMulti sources")
        
        # Use global DASH_DATA_PATH or fall back to legacy path
        import os
        data_root = os.environ.get('DASH_DATA_PATH', Path(__file__).parent.parent / "data")
        species_data_dir = Path(data_root) / self.species
        
        # Set BulkMulti dataset path
        self.data_dir = species_data_dir / "BulkMulti" / dataset
        self.atac_signals_dir = self.data_dir / "atac_signals"
        
        # Update database paths  
        self.bulkmulti_db_path = self.data_dir / "BulkMulti.db"
        
        # Reload sample configuration for new dataset
        self._load_sample_configuration()
        self._generate_sample_colors()
        
    def get_available_samples(self) -> List[str]:
        """Get available samples from BigWig files dynamically"""
        if not self.atac_signals_dir or not self.atac_signals_dir.exists():
            self.logger.warning(f"ATAC signals directory not found: {self.atac_signals_dir}")
            return self.sample_order  # Fallback to default
        
        try:
            # Get all BigWig files in the atac_signals directory
            bigwig_files = list(self.atac_signals_dir.glob("*.bw"))
            if not bigwig_files:
                self.logger.warning(f"No BigWig files found in {self.atac_signals_dir}")
                return self.sample_order  # Fallback to default
            
            # Extract sample names from file names
            sample_names = [bw_file.stem for bw_file in bigwig_files]  # Remove .bw extension
            
            # Sort samples to ensure consistent order
            sample_names.sort()
            
            return sample_names
            
        except Exception as e:
            self.logger.warning(f"Could not get samples from BigWig files: {e}")
            return self.sample_order  # Fallback to default
    
    def _load_sample_configuration(self):
        """Load sample configuration from TSV file or fallback to auto-detection"""
        if not self.data_dir:
            return
            
        config_file = self.data_dir / "sample_config.tsv"
        if config_file.exists():
            self.sample_config = self._parse_sample_config_tsv(config_file)
            self.sample_order = self.sample_config.get('order', [])
            self.logger.info(f"✅ Loaded sample configuration from {config_file.name}: {len(self.sample_order)} samples")
        else:
            # Fallback to auto-detection
            self.auto_detect_samples()
    
    def _parse_sample_config_tsv(self, config_file: Path) -> Dict:
        """Parse TSV format sample configuration file"""
        samples = []
        sample_colors = {}
        sample_display_names = {}
        sample_descriptions = {}
        
        try:
            with open(config_file, 'r', encoding='utf-8') as f:
                # Skip comment lines
                lines = [line for line in f if not line.strip().startswith('#')]
                
            # Parse TSV using csv.DictReader
            from io import StringIO
            csv_data = StringIO(''.join(lines))
            reader = csv.DictReader(csv_data, delimiter='\t')
            
            for row in reader:
                sample_name = row['sample_name'].strip()
                samples.append(sample_name)
                
                # Optional fields
                if row.get('display_name', '').strip():
                    sample_display_names[sample_name] = row['display_name'].strip()
                if row.get('color', '').strip():
                    sample_colors[sample_name] = row['color'].strip()
                if row.get('description', '').strip():
                    sample_descriptions[sample_name] = row['description'].strip()
            
            return {
                'order': samples,
                'colors': sample_colors,
                'display_names': sample_display_names,
                'descriptions': sample_descriptions
            }
            
        except Exception as e:
            self.logger.error(f"Error parsing sample config file {config_file}: {e}")
            return {}
    
    def auto_detect_samples(self):
        """Auto-detect samples from BigWig files and update sample_order"""
        detected_samples = self.get_available_samples()
        if detected_samples:
            self.sample_order = detected_samples
            self.sample_config = {'order': detected_samples, 'colors': {}, 'display_names': {}, 'descriptions': {}}
            self.logger.info(f"Auto-detected {len(detected_samples)} samples from BigWig files")
            # Regenerate colors for new samples
            self._generate_sample_colors()
    
    def _generate_sample_colors(self):
        """Generate colors for all samples dynamically"""
        # Initialize sample colors dictionary
        self.sample_colors = {}
        
        # Get colors from configuration if available
        config_colors = self.sample_config.get('colors', {})
        
        # Generate colors using Plotly color schemes for auto-assignment
        color_schemes = [
            px.colors.qualitative.Set1,
            px.colors.qualitative.Set2, 
            px.colors.qualitative.Set3,
            px.colors.qualitative.Pastel1,
            px.colors.qualitative.Dark2
        ]
        
        # Flatten all color schemes
        all_colors = []
        for scheme in color_schemes:
            all_colors.extend(scheme)
        
        # Assign colors to samples
        for i, sample in enumerate(self.sample_order):
            if sample in config_colors:
                # Use color from configuration file
                self.sample_colors[sample] = config_colors[sample]
            else:
                # Auto-assign color from Plotly palette (cycling)
                color_index = i % len(all_colors)
                self.sample_colors[sample] = all_colors[color_index]

    def _normalize_gene_list(self, gene_names: Optional[Sequence[str]]) -> List[str]:
        """Ensure provided gene names are valid strings."""
        if gene_names is None:
            return []
        if isinstance(gene_names, str):
            cleaned = gene_names.strip()
            return [cleaned] if cleaned else []
        if not isinstance(gene_names, Sequence):
            raise TypeError("gene_names must be a sequence of strings or None")
        normalized: List[str] = []
        for gene in gene_names:
            if gene is None:
                continue
            if not isinstance(gene, str):
                raise TypeError(f"Gene name must be a string, got {type(gene)!r}")
            stripped = gene.strip()
            if stripped:
                normalized.append(stripped)
        return normalized

    def _coerce_coordinate(self, value: Any, label: str) -> int:
        """Convert coordinate-like input to integer and validate."""
        if isinstance(value, bool):
            raise TypeError(f"{label} must not be boolean")
        try:
            coordinate = int(value)
        except (TypeError, ValueError) as exc:
            raise TypeError(f"{label} must be an integer") from exc
        if coordinate < 0:
            raise ValueError(f"{label} must be non-negative")
        return coordinate

    def _coerce_region_arguments(self, chromosome: str, start: Any, end: Any) -> Tuple[str, int, int]:
        """Validate chromosome/start/end inputs shared by region queries."""
        if not isinstance(chromosome, str):
            raise TypeError("chromosome must be a string")
        chrom = chromosome.strip()
        if not chrom:
            raise ValueError("chromosome cannot be empty")
        start_val = self._coerce_coordinate(start, "region_start")
        end_val = self._coerce_coordinate(end, "region_end")
        if start_val >= end_val:
            raise ValueError("region_start must be less than region_end")
        return chrom, start_val, end_val

    def _normalize_value_threshold(self, value_threshold: float) -> float:
        """Ensure linkage filters stay within valid ranges."""
        if isinstance(value_threshold, bool):
            raise TypeError("value_threshold cannot be boolean")
        try:
            threshold = float(value_threshold)
        except (TypeError, ValueError) as exc:
            raise TypeError("value_threshold must be numeric") from exc
        if threshold < 0:
            raise ValueError("value_threshold must be >= 0")
        return threshold
    
    def get_sample_color(self, sample_name: str) -> str:
        """Get color for a sample, with fallback generation"""
        if sample_name in self.sample_colors:
            return self.sample_colors[sample_name]
        
        # Fallback: generate color based on hash
        import hashlib
        hash_obj = hashlib.md5(sample_name.encode())
        hash_hex = hash_obj.hexdigest()
        # Convert first 6 characters to color
        color = f"#{hash_hex[:6]}"
        self.sample_colors[sample_name] = color
        return color
    
    def get_sample_display_name(self, sample_name: str) -> str:
        """Get display name for a sample, fallback to sample_name if not configured"""
        display_names = self.sample_config.get('display_names', {})
        return display_names.get(sample_name, sample_name)
    
    def get_sample_description(self, sample_name: str) -> str:
        """Get description for a sample"""
        descriptions = self.sample_config.get('descriptions', {})
        return descriptions.get(sample_name, '')
        
    def get_atac_data(self, chromosome: str, start: int, end: int) -> pd.DataFrame:
        """Get ATAC-seq data for genomic region from BigWig files"""
        if not self.atac_signals_dir or not self.atac_signals_dir.exists():
            self.logger.error(f"ATAC signals directory not found: {self.atac_signals_dir}")
            return pd.DataFrame()
        
        chromosome, start, end = self._coerce_region_arguments(chromosome, start, end)
        
        try:
            # Get all BigWig files in the atac_signals directory
            bigwig_files = list(self.atac_signals_dir.glob("*.bw"))
            if not bigwig_files:
                self.logger.error(f"No BigWig files found in {self.atac_signals_dir}")
                return pd.DataFrame()
            
            # Detect chromosome format from the first BigWig file
            test_bw = pyBigWig.open(str(bigwig_files[0]))
            available_chroms = list(test_bw.chroms().keys())
            test_bw.close()
            
            # Determine the correct chromosome format to use
            target_chromosome = None
            original_chrom = chromosome
            
            # Try different chromosome formats
            possible_formats = [
                chromosome,  # Use as-is
                f"chr{chromosome}" if not chromosome.startswith('chr') else chromosome[3:],  # Toggle chr prefix
                chromosome.replace('chr', '') if chromosome.startswith('chr') else f"chr{chromosome}"  # Alternative toggle
            ]
            
            for chrom_format in possible_formats:
                if chrom_format in available_chroms:
                    target_chromosome = chrom_format
                    break
            
            if not target_chromosome:
                self.logger.error(f"Chromosome {original_chrom} not found in any format. Available: {available_chroms[:5]}...")
                return pd.DataFrame()
            
            self.logger.info(f"Using chromosome format '{target_chromosome}' for region {original_chrom}")
            
            # Create data structure for all samples
            atac_data = []
            
            for bw_file in bigwig_files:
                sample_name = bw_file.stem  # Remove .bw extension
                
                try:
                    # Open BigWig file
                    bw = pyBigWig.open(str(bw_file))
                    
                    # Check if chromosome exists in this BigWig file
                    if target_chromosome not in bw.chroms():
                        self.logger.warning(f"Chromosome {target_chromosome} not found in {bw_file}")
                        bw.close()
                        continue
                    
                    # Get values for the region (bin size = 100bp for visualization)
                    bin_size = 100
                    total_length = max(end - start, bin_size)
                    n_bins = max(1, int(np.ceil(total_length / bin_size)))
                    stats = bw.stats(
                        target_chromosome,
                        start,
                        end,
                        nBins=n_bins,
                        type='mean',
                        exact=False
                    ) or []
                    
                    for idx in range(n_bins):
                        bin_start = start + idx * bin_size
                        bin_end = min(bin_start + bin_size, end)
                        value = stats[idx] if idx < len(stats) else None
                        
                        if value is None or pd.isna(value):
                            avg_signal = 0.0
                        else:
                            avg_signal = float(value)
                        
                        atac_data.append({
                            'chr': target_chromosome,
                            'start': bin_start,
                            'end': bin_end,
                            'sample': sample_name,
                            'signal': avg_signal
                        })
                    
                    bw.close()
                    
                except Exception as e:
                    self.logger.error(f"Error reading BigWig file {bw_file}: {e}")
                    continue
            
            if not atac_data:
                self.logger.warning(f"No ATAC data found for region {chromosome}:{start}-{end}")
                return pd.DataFrame()
            
            # Convert to DataFrame and pivot for compatibility with existing code
            df = pd.DataFrame(atac_data)
            
            # Pivot to get samples as columns (similar to old SQLite format)
            pivot_df = df.pivot_table(
                index=['chr', 'start', 'end'], 
                columns='sample', 
                values='signal', 
                fill_value=0.0
            ).reset_index()
            
            # Flatten column names
            pivot_df.columns = [col if col in ['chr', 'start', 'end'] else col for col in pivot_df.columns]
            
            return pivot_df
            
        except Exception as e:
            self.logger.error(f"Error getting ATAC data from BigWig files: {e}")
            return pd.DataFrame()
    
    def get_rna_data(self, gene_names: Optional[Sequence[str]]) -> pd.DataFrame:
        """Get RNA expression data for genes from BulkMulti.db"""
        if not self.bulkmulti_db_path or not self.bulkmulti_db_path.exists():
            self.logger.error(f"BulkMulti database not found: {self.bulkmulti_db_path}")
            return pd.DataFrame()
        
        normalized_genes = self._normalize_gene_list(gene_names)
        if not normalized_genes:
            return pd.DataFrame()
        
        try:
            with self._connect_ro(self.bulkmulti_db_path) as conn:
                # Check available tables
                tables_df = pd.read_sql_query("SELECT name FROM sqlite_master WHERE type='table'", conn)
                available_tables = tables_df['name'].tolist()
                # Try different possible table names for RNA data
                possible_tables = ['rna_expression', 'Expression_TPM', 'rna_tpm', 'gene_expression']
                rna_table = None

                for table in possible_tables:
                    if table in available_tables:
                        rna_table = table
                        break

                if not rna_table:
                    self.logger.error(f"No RNA tables found in BulkMulti.db. Available tables: {available_tables}")
                    return pd.DataFrame()

                # Create placeholders for gene names
                placeholders = ','.join(['?'] * len(normalized_genes))

                # Try different possible column names for gene names
                possible_name_cols = ['name', 'gene_name', 'ens_id', 'gene_id']

                for name_col in possible_name_cols:
                    try:
                        query = f"SELECT * FROM {rna_table} WHERE {name_col} IN ({placeholders})"
                        df = pd.read_sql_query(query, conn, params=normalized_genes)
                        if not df.empty:
                            self.logger.info(f"✅ RNA data loaded from {rna_table} table using {name_col} column")
                            return df
                    except Exception as exc:
                        self.logger.debug(
                            "Failed RNA query for table=%s column=%s genes=%s: %s",
                            rna_table,
                            name_col,
                            normalized_genes,
                            exc,
                        )

                self.logger.warning(
                    "Found RNA table %s but could not match gene names %s using columns %s",
                    rna_table,
                    normalized_genes,
                    possible_name_cols,
                )
                return pd.DataFrame()
            
        except Exception as e:
            self.logger.error(f"Error getting RNA data from BulkMulti.db: {e}")
            return pd.DataFrame()
    
    def create_atac_tracks(self, atac_data: pd.DataFrame, region_start: int, region_end: int) -> List[go.Scatter]:
        """Create ATAC-seq tracks similar to the reference image"""
        traces = []
        
        if atac_data.empty:
            return traces
        
        # Create one trace for each sample
        for i, sample in enumerate(self.sample_order):
            if sample in atac_data.columns:
                # Sort data by genomic position
                sorted_data = atac_data.sort_values('start')
                
                # Create step-like visualization for each bin
                x_coords = []
                y_coords = []
                
                for _, row in sorted_data.iterrows():
                    # Create rectangular area for each bin
                    # Ensure numeric values and handle NaN/infinity
                    start_pos = int(row['start'])
                    end_pos = int(row['end'])
                    signal_val = float(row[sample])
                    
                    # Skip if signal is NaN or infinite
                    if pd.isna(signal_val) or not np.isfinite(signal_val):
                        signal_val = 0.0
                    
                    x_coords.extend([start_pos, start_pos, end_pos, end_pos])
                    y_coords.extend([0.0, signal_val, signal_val, 0.0])
                
                # Only create trace if we have data
                if len(x_coords) > 0:
                    # Add transparency to fill color using dynamic color generation
                    base_color = self.get_sample_color(sample)
                    if base_color.startswith('#'):
                        # Convert hex to rgba with transparency
                        hex_color = base_color.lstrip('#')
                        rgb = tuple(int(hex_color[i:i+2], 16) for i in (0, 2, 4))
                        fill_color = f'rgba({rgb[0]}, {rgb[1]}, {rgb[2]}, 0.7)'
                    else:
                        fill_color = base_color
                    
                    trace = go.Scatter(
                        x=x_coords,
                        y=y_coords,
                        fill='tozeroy',
                        fillcolor=fill_color,
                        line=dict(color=base_color, width=0.5),
                        name=self.get_sample_display_name(sample),
                        showlegend=False,
                        hovertemplate=f"{self.get_sample_display_name(sample)}<br>Position: %{{x:,}}<br>Signal: %{{y:.3f}}<extra></extra>"
                    )
                    traces.append(trace)
        
        return traces
    
    def create_expression_plot(self, rna_data: pd.DataFrame, gene_name: str) -> go.Bar:
        """Create expression bar plot for a gene"""
        if rna_data.empty:
            return None
        
        gene_row = rna_data[rna_data['name'] == gene_name]
        if gene_row.empty:
            return None
        
        # Extract expression values for each sample in reverse order to match ATAC tracks
        expression_values = []
        sample_names = []
        
        for sample in self.sample_order:
            if sample in gene_row.columns:
                expression_values.append(gene_row[sample].iloc[0])
                sample_names.append(self.get_sample_display_name(sample))
        
        # Reverse the order to match ATAC tracks (top to bottom)
        expression_values.reverse()
        sample_names.reverse()
        
        # Create horizontal bar plot with proper spacing
        trace = go.Bar(
            y=sample_names,
            x=expression_values,
            orientation='h',
            marker_color=[self.get_sample_color(sample) for sample in self.sample_order if sample in gene_row.columns][::-1],
            name=gene_name,
            showlegend=False,
            width=0.8,  # Reduce bar width for better spacing
            hovertemplate=f"{gene_name}<br>%{{y}}<br>Expression: %{{x:.2f}} TPM<extra></extra>"
        )
        
        return trace
    
    def create_linkage_arcs(self, linkage_data: pd.DataFrame, region_start: int, region_end: int) -> List[go.Scatter]:
        """Create linkage arcs similar to the reference image"""
        traces = []
        
        if linkage_data.empty:
            return traces
        
        for _, linkage in linkage_data.iterrows():
            # Calculate arc parameters
            x_start = linkage['Start']
            x_end = linkage['End']
            x_center = (x_start + x_end) / 2
            arc_height = (x_end - x_start) / 10  # Adjust arc height based on distance
            
            # Create arc points
            t = np.linspace(0, np.pi, 50)
            x_arc = x_center + (x_end - x_start) / 2 * np.cos(t)
            y_arc = -arc_height * np.sin(t)
            
            trace = go.Scatter(
                x=x_arc,
                y=y_arc,
                mode='lines',
                line=dict(
                    color=px.colors.sequential.Viridis[int(linkage['Value'] * 10)],
                    width=2
                ),
                showlegend=False,
                hovertemplate=f"Linkage<br>Value: {linkage['Value']:.3f}<br>Distance: {linkage['Dist']:,} bp<extra></extra>"
            )
            traces.append(trace)
        
        return traces
    
    def create_gene_track(self, gene_data: pd.DataFrame, region_start: int, region_end: int) -> List[go.Scatter]:
        """Create gene annotation track with proper gene structure"""
        traces = []
        
        if gene_data.empty:
            return traces
        
        y_level = 0.5
        for _, gene in gene_data.iterrows():
            # Gene body as rectangle
            gene_trace = go.Scatter(
                x=[gene['Start'], gene['End'], gene['End'], gene['Start'], gene['Start']],
                y=[y_level-0.2, y_level-0.2, y_level+0.2, y_level+0.2, y_level-0.2],
                mode='lines',
                fill='toself',
                fillcolor='red' if gene['Strand'] == '+' else 'blue',
                line=dict(
                    color='red' if gene['Strand'] == '+' else 'blue',
                    width=2
                ),
                name=gene['Symbol'],
                showlegend=False,
                hovertemplate=f"{gene['Symbol']}<br>Position: {gene['Start']:,}-{gene['End']:,}<br>Strand: {gene['Strand']}<extra></extra>"
            )
            traces.append(gene_trace)
            
            # Gene label - positioned based on strand
            label_x = gene['Start'] if gene['Strand'] == '+' else gene['End']
            label_trace = go.Scatter(
                x=[label_x],
                y=[y_level + 0.7],
                mode='text',
                text=[gene['Symbol']],
                textposition='middle center',
                textfont=dict(size=10, color='black'),
                showlegend=False,
                hoverinfo='skip'
            )
            traces.append(label_trace)
            
            y_level += 1  # Stack genes vertically if multiple
        
        return traces
    
    
    def create_landscape_view_plot(
        self,
        gene_names: Optional[Sequence[str]],
        chromosome: str,
        region_start: int,
        region_end: int,
        bulk_multi_api=None,
        value_threshold: float = 0.2,
        render_config: Optional[Dict[str, Any]] = None
    ) -> go.Figure:
        """Create the main landscape view plot matching the reference image"""
        try:
            normalized_genes = self._normalize_gene_list(gene_names)
            chromosome, region_start, region_end = self._coerce_region_arguments(chromosome, region_start, region_end)
            linkage_threshold = self._normalize_value_threshold(value_threshold)
            region_start_display = f"{region_start:,}"
            region_end_display = f"{region_end:,}"
            # Get all required data
            atac_data = self.get_atac_data(chromosome, region_start, region_end)
            rna_data = self.get_rna_data(normalized_genes)
            settings = self._get_render_settings(render_config)
            font_family = settings['font_family']
            panel_title_font_size = settings['panel_title_font_size']
            track_label_font_size = settings['track_label_font_size']
            gene_label_font_size = settings['gene_label_font_size']
            atac_signal_ratio = settings['atac_height_ratio']
            show_els_lines = settings['show_els_guide_lines']
            atac_track_height = settings['atac_track_height']
            gene_track_height = settings['gene_track_height']
            els_track_height = settings['els_track_height']
            linkage_track_height = settings['linkage_track_height']
            panel_title_height = settings['panel_title_height']
            
            # Debug ATAC data retrieval
            # Basic data validation
            if atac_data.empty:
                self.logger.warning("⚠️ No ATAC data found for this region")
            
            # Auto-detect samples if not explicitly set
            if not self.sample_order:
                self.auto_detect_samples()
                self.logger.debug("📊 Auto-detected %d ATAC samples", len(self.sample_order))
            
            # Initialize gene annotations list
            gene_annotations = []
            
            # Create improved subplot structure: ATAC + genes + spacing + ELS + linkages
            n_atac_tracks = len(self.sample_order)
            spacing_height = 20  # Keep spacing between genes and ELS for clarity
            total_rows = n_atac_tracks + 4  # ATAC tracks + genes + spacing + ELS + linkages
            
            # Calculate heights based on fixed pixel heights
            total_track_height = (n_atac_tracks * atac_track_height + 
                                gene_track_height + 
                                spacing_height + 
                                els_track_height + 
                                linkage_track_height)
            
            # Calculate total figure height including margins and titles
            total_figure_height = total_track_height + panel_title_height + 100  # 100px for margins
            
            # Convert pixel heights to proportions
            atac_height_ratio = atac_track_height / total_track_height
            gene_height_ratio = gene_track_height / total_track_height
            spacing_height_ratio = spacing_height / total_track_height
            els_height_ratio = els_track_height / total_track_height
            linkage_height_ratio = linkage_track_height / total_track_height
            
            row_heights = [atac_height_ratio] * n_atac_tracks + [gene_height_ratio, spacing_height_ratio, els_height_ratio, linkage_height_ratio]
            
            
            fig = make_subplots(
                rows=total_rows,
                cols=2,
                shared_xaxes=True,
                vertical_spacing=0,  # Reduced for more compact layout
                horizontal_spacing=0,  # Reduced horizontal spacing between panels
                column_widths=[0.88, 0.12],  # Make expression panel much narrower
                row_heights=row_heights,
                specs=[[{"type": "scatter"}, {"type": "bar", "rowspan": n_atac_tracks}] if i == 0 
                       else [{"type": "scatter"}, None] 
                       for i in range(total_rows)]
            )
            
            # Add ATAC-seq tracks
            atac_traces = self.create_atac_tracks(atac_data, region_start, region_end)
            for i, trace in enumerate(atac_traces):
                fig.add_trace(trace, row=i+1, col=1)
            
            # Add expression plot (right panel) - spans all ATAC track rows
            if not rna_data.empty and normalized_genes:
                gene_name = normalized_genes[0]  # Use first gene
                exp_trace = self.create_expression_plot(rna_data, gene_name)
                if exp_trace:
                    fig.add_trace(exp_trace, row=1, col=2)
            
            # Get ELS positions from peaks table and linkages separately, then match them
            els_y_pos = 0.5
            feature_positions = []
            region_linkages = pd.DataFrame()  # Store filtered linkages for later use
            
            # Get real ELS (peaks) and linkage data from bulk_multi_api
            if bulk_multi_api:
                try:
                    # Step 1: Get actual ELS positions from peaks table in the region
                    import sqlite3
                    peaks_query = """
                        SELECT Chromosome, Start, End FROM peaks 
                        WHERE Chromosome = ? 
                        AND Start >= ? AND Start <= ?
                        ORDER BY Start
                    """
                    with self._connect_ro(bulk_multi_api.db_path) as conn:
                        peaks_df = pd.read_sql_query(
                            peaks_query, conn, 
                            params=[str(chromosome), int(region_start), int(region_end)]
                        )
                    self.logger.debug("📍 Found %d ELS peaks in %s:%s-%s", len(peaks_df), chromosome, region_start_display, region_end_display)
                    
                    if not peaks_df.empty:
                        # Use peak centers as ELS positions
                        feature_positions = sorted([
                            int((row['Start'] + row['End']) / 2) 
                            for _, row in peaks_df.iterrows()
                        ])
                    
                    # Step 2: Get linkages for current genes
                    linkages_df = bulk_multi_api.get_gene_linkages(
                        normalized_genes, linkage_threshold
                    )
                    if not linkages_df.empty:
                        # Filter linkages: require gene (End) in region
                        region_linkages = linkages_df[
                            (linkages_df['Chromosome'] == chromosome) &
                            (linkages_df['End'] >= region_start) & 
                            (linkages_df['End'] <= region_end)
                        ]
                        
                        # Step 3: Further filter linkages to only show those connecting to actual ELS positions
                        if feature_positions and not region_linkages.empty:
                            # Find linkages that connect to actual peaks (within reasonable distance)
                            matched_linkages = []
                            for _, linkage in region_linkages.iterrows():
                                linkage_start = linkage['Start']
                                # Check if this linkage connects to any actual peak (within 1kb tolerance)
                                for peak_pos in feature_positions:
                                    if abs(linkage_start - peak_pos) <= 1000:
                                        matched_linkages.append(linkage)
                                        break
                            
                            if matched_linkages:
                                region_linkages = pd.DataFrame(matched_linkages)
                                self.logger.debug("🔗 Found %d gene-ELS linkages in region", len(region_linkages))
                            else:
                                region_linkages = pd.DataFrame()
                        
                except Exception as e:
                    self.logger.warning(f"⚠️ Could not get ELS/linkage data: {e}")
            
            # If no real ELS positions found, continue without ELS features
            if not feature_positions:
                pass  # Continue without ELS features
            
            # Add ELS guide lines to ATAC tracks if enabled
            if show_els_lines and feature_positions:
                for pos in feature_positions:
                    for track_idx in range(n_atac_tracks):
                        fig.add_vline(
                            x=pos,
                            line=dict(color="green", width=1, dash="dash"),
                            row=track_idx + 1, col=1
                        )
            
            # Add ELS feature ticks in ELS track for each regulatory element
            for pos in feature_positions:
                fig.add_trace(
                    go.Scatter(
                        x=[pos, pos],
                        y=[0, 1],
                        mode='lines',
                        line=dict(color='green', width=3),
                        showlegend=False,
                        hovertemplate=f"Regulatory Element<br>Position: {pos:,}<extra></extra>"
                    ),
                    row=n_atac_tracks + 3, col=1  # ELS track
                )
            
            # Add linkage arcs using the SAME filtered data from above
            if not region_linkages.empty:
                try:
                    
                    # Create arcs for each linkage with improved visibility
                    for i, (_, linkage) in enumerate(region_linkages.iterrows()):
                        if i > 15:  # Limit to prevent overcrowding
                            break
                            
                        x_start = linkage['Start']
                        x_end = linkage['End']
                        x_center = (x_start + x_end) / 2
                        arc_height = min((x_end - x_start) / 12000, 1.8)  # Improved arc height
                        
                        # Create arc points with fewer points for better performance
                        t = np.linspace(0, np.pi, 20)  # Reduced from 60 to 20 points for speed
                        x_arc = x_center + (x_end - x_start) / 2 * np.cos(t)
                        y_arc = -arc_height * np.sin(t)  # Negative for downward arcs (inverted)
                        
                        # Enhanced color scheme based on linkage strength
                        color_intensity = min(linkage['Value'], 1.0)
                        if color_intensity > 0.7:
                            arc_color = f'rgba(255, 0, 0, {0.6 + color_intensity * 0.4})'  # Red for strong
                        elif color_intensity > 0.4:
                            arc_color = f'rgba(255, 165, 0, {0.5 + color_intensity * 0.5})'  # Orange for medium
                        else:
                            arc_color = f'rgba(100, 149, 237, {0.4 + color_intensity * 0.6})'  # Blue for weak
                        
                        fig.add_trace(
                            go.Scatter(
                                x=x_arc,
                                y=y_arc,
                                mode='lines',
                                line=dict(color=arc_color, width=2.5),
                                showlegend=False,
                                hovertemplate=f"Linkage Strength: {linkage['Value']:.3f}<br>Start: {x_start:,} bp<br>End: {x_end:,} bp<br>Distance: {abs(x_end-x_start):,} bp<extra></extra>"
                            ),
                            row=n_atac_tracks + 4, col=1  # Linkages are now at bottom
                        )
                except Exception as e:
                    self.logger.warning(f"Could not add linkage arcs: {e}")
            
            # Add gene track with real gene data
            if bulk_multi_api:
                try:
                    # Get all genes in the region, not just target genes
                    genes_in_region = bulk_multi_api.search_genes_by_region(
                        chromosome, region_start, region_end
                    )
                    
                    # Also include target genes if they're not in the region
                    all_gene_names = list(set(genes_in_region + normalized_genes))
                    
                    if all_gene_names:
                        self.logger.debug("🧬 Found %d genes in region %s:%s-%s", len(genes_in_region), chromosome, region_start_display, region_end_display)
                        if normalized_genes:
                            self.logger.debug("🎯 Target genes: %s", normalized_genes)
                        self.logger.debug("📊 Total genes to display: %d", len(all_gene_names))
                        gene_info = bulk_multi_api.get_gene_info(all_gene_names)
                        if not gene_info.empty:
                            self.logger.debug("📊 Displaying gene structures for %d genes", len(gene_info))
                            # Batch fetch all exon data for performance
                            all_exons = bulk_multi_api.get_exon_data(
                                all_gene_names, chromosome, region_start, region_end
                            )
                            
                            # Sort genes by start position for better layout
                            gene_info_sorted = gene_info.sort_values('Start')
                            y_level = 0.5
                            for _, gene in gene_info_sorted.iterrows():
                                # Gene body as rectangle
                                gene_start = gene['Start']
                                gene_end = gene['End']
                                
                                # Only show genes that overlap with the region
                                if gene_end >= region_start and gene_start <= region_end:
                                    strand = gene.get('Strand', '+')
                                    gene_color = 'red' if strand == '+' else 'blue'
                                    
                                    # Gene structure: intron line + exon boxes
                                    visible_start = max(gene_start, region_start)
                                    visible_end = min(gene_end, region_end)
                                    
                                    # Draw intron line without hover
                                    intron_trace = go.Scatter(
                                        x=[visible_start, visible_end],
                                        y=[y_level, y_level],
                                        mode='lines',
                                        line=dict(color=gene_color, width=2),
                                        showlegend=False,
                                        hoverinfo='skip'  # No hover for intron
                                    )
                                    fig.add_trace(intron_trace, row=n_atac_tracks + 1, col=1)
                                    
                                    # Get real exon data for this gene from pre-fetched data
                                    gene_symbol = gene['Symbol']
                                    try:
                                        real_exons = all_exons[all_exons['Gene_Symbol'] == gene_symbol]
                                        
                                        if not real_exons.empty:
                                            # Draw real exons as shapes
                                            exon_height = 0.4
                                            for _, exon in real_exons.iterrows():
                                                exon_start = exon['Start']
                                                exon_end = exon['End']
                                                
                                                # Only draw if exon is visible in region
                                                if exon_end >= region_start and exon_start <= region_end:
                                                    visible_exon_start = max(exon_start, region_start)
                                                    visible_exon_end = min(exon_end, region_end)
                                                    
                                                    fig.add_shape(
                                                        type="rect",
                                                        x0=visible_exon_start,
                                                        x1=visible_exon_end,
                                                        y0=y_level - exon_height,
                                                        y1=y_level + exon_height,
                                                        fillcolor=gene_color,
                                                        line=dict(color=gene_color, width=1),
                                                        row=n_atac_tracks + 1, col=1
                                                    )
                                            
                                        else:
                                            # Fallback to simple gene rectangle if no exons found
                                            pass
                                            
                                    except Exception as e:
                                        # Fallback - no exon display, just intron line
                                        pass
                                    
                                    # Add invisible trace for gene hover over entire gene region
                                    gene_hover_trace = go.Scatter(
                                        x=[visible_start, visible_end],
                                        y=[y_level, y_level],
                                        mode='lines',
                                        line=dict(color='rgba(0,0,0,0)', width=20),  # Invisible but wide for easy hovering
                                        showlegend=False,
                                        hovertemplate=f"{gene['Symbol']}<br>Start: {gene_start:,} bp<br>End: {gene_end:,} bp<br>Strand: {strand}<br>Length: {gene_end-gene_start:,} bp<extra></extra>"
                                    )
                                    fig.add_trace(gene_hover_trace, row=n_atac_tracks + 1, col=1)
                                    
                                    # Strand arrows disabled for cleaner visualization
                                    # (Gene strand info is still shown by gene name color: red=+, blue=-)
                                    
                                    # Gene label with improved positioning to avoid overlap
                                    label_x = (visible_start + visible_end) / 2
                                    
                                    # Add gene label directly to the gene track using scatter plot
                                    # This keeps the label within the track coordinate system
                                    gene_label_trace = go.Scatter(
                                        x=[label_x],
                                        y=[y_level + 1.5],  # Slightly above the gene center
                                        mode='text',
                                        text=[gene['Symbol']],
                                        textposition='middle center',
                                        textfont=dict(size=gene_label_font_size, color='blue' if gene.get('Strand', '+') == '-' else 'red', family=font_family),
                                        showlegend=False,
                                        hoverinfo='skip'
                                    )
                                    fig.add_trace(gene_label_trace, row=n_atac_tracks + 1, col=1)
                                    
                                    y_level += 1.0  # Increased spacing to avoid overlap
                except Exception as e:
                    self.logger.warning(f"⚠️ Could not add gene track: {e}")
                    # Fallback to simple label within track
                    if normalized_genes:
                        gene_name = normalized_genes[0]
                        fallback_trace = go.Scatter(
                            x=[(region_start + region_end) / 2],
                            y=[0.5],
                            mode='text',
                            text=[f"Gene: {gene_name}"],
                            textposition='middle center',
                            textfont=dict(size=gene_label_font_size, color="red", family=font_family),
                            showlegend=False,
                            hoverinfo='skip'
                        )
                        fig.add_trace(fallback_trace, row=n_atac_tracks + 1, col=1)
            
            # Update x-axis for the bottom subplot only
            fig.update_xaxes(
                range=[region_start, region_end],
                title_text="",  # Remove x-axis title
                showline=True,  # Show x-axis line
                linewidth=2,  # Thicker line
                linecolor='black',
                showticklabels=True,
                tickfont=dict(size=10),
                tickformat=',d',  # Format numbers with commas
                ticks="outside",  # Show tick marks outside the plot
                tickwidth=1,  # Tick mark width
                ticklen=5,  # Tick mark length
                tickcolor='black',  # Tick mark color
                row=total_rows, col=1
            )
            
            # Initialize annotations list for all labels
            annotations = []
            
            # Calculate cumulative heights for proper positioning based on actual row heights
            cumulative_height = 0
            track_positions = []
            for i, height in enumerate(row_heights):
                # Position label at center of each track
                y_position = 1.0 - cumulative_height - (height / 2)
                track_positions.append(y_position)
                cumulative_height += height
            
            # Update y-axis labels and scaling for each ATAC track using actual heights
            for i, sample in enumerate(self.sample_order):
                max_signal = atac_data[sample].max() if not atac_data.empty and sample in atac_data.columns else 1
                y_position = track_positions[i]
                
                annotations.append(
                    dict(
                        x=-0.03,  # Left side of plot area
                        y=y_position,  # Position based on actual row height
                        xref='paper',
                        yref='paper',
                        text=self.get_sample_display_name(sample),
                        showarrow=False,
                        font=dict(size=track_label_font_size, color='black', family=font_family),
                        xanchor='right',
                        yanchor='middle',
                        textangle=0  # Horizontal text
                    )
                )
                
                # Calculate y-axis range using the height ratio parameter
                if max_signal > 0:
                    # Find overall max signal across all samples for consistent scaling
                    overall_max = max([atac_data[s].max() for s in self.sample_order if s in atac_data.columns and not atac_data[s].empty])
                    display_max = overall_max * atac_signal_ratio
                else:
                    display_max = 1
                
                fig.update_yaxes(
                    row=i+1, 
                    col=1, 
                    showticklabels=False,
                    range=[0, display_max],  # Use consistent scaled range for all tracks
                    fixedrange=True,
                    showgrid=False
                )
                
                # Add x-axis line for each ATAC track
                fig.update_xaxes(
                    row=i+1,
                    col=1,
                    showline=True,
                    linewidth=0.5,
                    linecolor='lightgray'
                )
            
            # Note: ELS guide lines are handled earlier in the function based on show_els_guide_lines config
            
            # Add horizontal labels for other tracks using actual positions
            other_tracks = ["Genes", "ELS", "Linkages"]
            track_indices = [n_atac_tracks, n_atac_tracks + 2, n_atac_tracks + 3]  # Skip spacing row
            for j, track_name in enumerate(other_tracks):
                track_row_index = track_indices[j]
                y_position = track_positions[track_row_index]
                
                annotations.append(
                    dict(
                        x=-0.03,  # Left side of plot area
                        y=y_position,  # Position based on actual row height
                        xref='paper',
                        yref='paper',
                        text=track_name,
                        showarrow=False,
                        font=dict(size=track_label_font_size, color='black', family=font_family),
                        xanchor='right',
                        yanchor='middle',
                        textangle=0  # Horizontal text
                    )
                )
            
            # Set y-axis ranges without titles (genes, spacing, ELS, linkages)
            fig.update_yaxes(range=[0, 6], row=n_atac_tracks + 1, col=1, showticklabels=False)  # Genes track
            fig.update_yaxes(range=[0, 1], row=n_atac_tracks + 2, col=1, showticklabels=False, showgrid=False, visible=False)  # Spacing track (hidden)
            fig.update_yaxes(range=[0, 1], row=n_atac_tracks + 3, col=1, showticklabels=False)  # ELS track
            fig.update_yaxes(range=[-3, 1], row=n_atac_tracks + 4, col=1, showticklabels=False)  # Linkages track (negative for downward arcs)
            
            # Add x-axis line only for linkages track (bottom) - thicker with ticks
            fig.update_xaxes(
                row=n_atac_tracks + 4, col=1, 
                showline=True, 
                linewidth=2, 
                linecolor='black',
                showticklabels=True,
                tickfont=dict(size=10),
                tickformat=',d',  # Format numbers with commas
                ticks="outside",  # Show tick marks outside the plot
                tickwidth=1,  # Tick mark width
                ticklen=5,  # Tick mark length
                tickcolor='black'  # Tick mark color
            )
            
            # Set expression plot y-axis (no sample labels to make it more compact)
            fig.update_yaxes(row=1, col=2, showticklabels=False)
                    
            # Add panel titles and track labels as annotations
            # Note: annotations list was already started above for track labels
            
            # ATAC panel title (well above the plot area)
            annotations.append(
                dict(
                    x=0.44,  # Center of left pane
                    y=1.01,  # above the figure
                    xref='paper',
                    yref='paper',
                    text="Accessibility",
                    showarrow=False,
                    font=dict(size=panel_title_font_size, color='black', family=font_family),
                    xanchor='center',
                    yanchor='middle'
                )
            )
            
            # Expression panel title (well above the plot area)
            annotations.append(
                dict(
                    x=0.94,  # Center of right panel (adjusted for new width: 88% + 12% * 0.5)
                    y=1.01,   # Well above the figure
                    xref='paper',
                    yref='paper',
                    text="Expression (TPM)",
                    showarrow=False,
                    font=dict(size=panel_title_font_size, color='black', family=font_family),
                    xanchor='center',
                    yanchor='middle'
                )
            )
            
            # Combine track labels, panel titles, and gene annotations
            all_annotations = annotations + gene_annotations
            
            # Update layout for better visibility
            fig.update_layout(
                title=f"Genomic Browser View: {chromosome}:{int(region_start):,}-{int(region_end):,}",
                height=total_figure_height,  # Dynamic height based on track count and fixed heights
                showlegend=False,
                plot_bgcolor='white',
                margin=dict(l=160, r=30, t=120, b=60),  # Adjusted margins: reduced left, increased top
                annotations=all_annotations
            )
            
            return fig
            
        except Exception as e:
            self.logger.error(f"Error creating genomic browser plot: {e}")
            # Return empty figure
            return go.Figure()
