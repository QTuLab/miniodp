#!/usr/bin/env python3
"""
Landscape View CLI - command line utility built on BulkMulti data that
renders multi-omics landscape figures (PNG/HTML) for Peak2Gene analysis.
"""

import sys
import argparse
import logging
from pathlib import Path
from typing import Dict, List, Optional


SCRIPTS_DIR = Path(__file__).resolve().parent
DASH_DIR = SCRIPTS_DIR.parent
if str(DASH_DIR) not in sys.path:
    sys.path.insert(0, str(DASH_DIR))
if str(SCRIPTS_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPTS_DIR))

from logging_utils import configure_script_logger, create_print_proxy  # type: ignore

_SCRIPT_LOGGER = configure_script_logger(__name__)
print = create_print_proxy(_SCRIPT_LOGGER)


class LandscapeViewCLI:
    """Landscape View command line driver."""
    
    def __init__(self):
        self.logger = self._setup_logging()
        self.config = {}
        self.landscape_api = None
        self.bulk_api = None
    
    def _setup_logging(self, debug=False) -> logging.Logger:
        """Configure module logger."""
        level = logging.DEBUG if debug else logging.INFO
        logging.basicConfig(
            level=level,
            format='%(asctime)s - %(levelname)s - %(message)s',
            datefmt='%H:%M:%S'
        )
        return logging.getLogger(__name__)
    
    def load_config(self, config_path: str) -> bool:
        """Load TOML config."""
        try:
            try:
                import toml  # type: ignore
            except ImportError as exc:
                raise RuntimeError("Missing dependency: toml (install via `pip install toml`)") from exc

            config_file = Path(config_path)
            if not config_file.exists():
                self.logger.error(f"Config file not found: {config_path}")
                return False
            
            with open(config_file, 'r', encoding='utf-8') as f:
                self.config = toml.load(f)
            
            self.logger.info(f"✅ Loaded config file: {config_path}")
            return True
            
        except Exception as e:
            self.logger.error(f"Failed to load config file: {e}")
            return False
    
    def validate_config(self) -> bool:
        """Validate config content."""
        required_sections = ['general', 'target', 'thresholds', 'visualization', 'export']
        
        for section in required_sections:
            if section not in self.config:
                self.logger.error(f"Config file missing required section: {section}")
                return False
        
        # Validate required fields
        required_fields = {
            'general': ['species', 'output_dir', 'output_name', 'formats'],
            'target': ['gene_names', 'upstream_bp', 'downstream_bp'],
            'thresholds': ['value_threshold'],
            'export': ['png_width', 'png_height']
        }
        
        for section, fields in required_fields.items():
            for field in fields:
                if field not in self.config[section]:
                    self.logger.error(f"Section {section} missing field: {field}")
                    return False
        
        # Validate output directory
        output_dir = Path(self.config['general']['output_dir'])
        if not output_dir.exists():
            self.logger.error(f"Output directory does not exist: {output_dir}")
            return False
        
        # Validate gene list
        gene_names = self.config['target']['gene_names']
        if not gene_names or not isinstance(gene_names, list):
            self.logger.error("Gene name list must be a non-empty list")
            return False
        
        self.logger.info("✅ Config validated successfully")
        return True
    
    def initialize_apis(self) -> bool:
        """Initialize BulkMulti and visualization APIs."""
        try:
            from backend.bulk_multi_api import BulkMultiAPI
            from backend.landscape_view_api import LandscapeViewAPI

            species = self.config['general']['species']
            
            # Initialize BulkMulti data API
            self.bulk_api = BulkMultiAPI(species_key=species)
            
            # Build visualization configuration
            viz_config = self.config.get('visualization', {})
            custom_config = {
                'atac_height_ratio': viz_config.get('atac_height_ratio', 0.01),
                'atac_height': viz_config.get('atac_height', 60),
                'gene_height': viz_config.get('gene_height', 200),
                'els_height': viz_config.get('els_height', 20),
                'linkage_height': viz_config.get('linkage_height', 300),
                'font_family': viz_config.get('font_family', 'Arial'),
                'panel_title_size': viz_config.get('panel_title_size', 18),
                'track_label_size': viz_config.get('track_label_size', 18),
                'gene_label_size': viz_config.get('gene_label_size', 16),
                'sample_order': viz_config.get('sample_order', []),
                'show_els_guide_lines': viz_config.get('show_els_guide_lines', False)
            }
            
            # Initialize visualization API
            self.landscape_api = LandscapeViewAPI(species_key=species, custom_config=custom_config)
            
            # Auto-detect sample order when not supplied
            if not custom_config['sample_order']:
                self.landscape_api.auto_detect_samples()
            
            self.logger.info("✅ APIs initialized")
            return True
            
        except Exception as e:
            self.logger.error(f"Failed to initialize APIs: {e}")
            return False
    
    def calculate_region(self) -> tuple:
        """Derive visualization genomic region based on target genes."""
        gene_names = self.config['target']['gene_names']
        upstream_bp = self.config['target']['upstream_bp']
        downstream_bp = self.config['target']['downstream_bp']
        
        # Fetch gene coordinates
        gene_info = self.bulk_api.get_gene_info(gene_names)
        if gene_info.empty:
            raise ValueError(f"Gene metadata unavailable: {gene_names}")
        
        # Determine gene interval boundaries
        gene_region_start = min(gene_info['Start'])
        gene_region_end = max(gene_info['End'])
        
        # Expand region by upstream/downstream padding
        region_start = gene_region_start - upstream_bp
        region_end = gene_region_end + downstream_bp
        
        self.logger.info(f"Gene interval: {gene_region_start:,}-{gene_region_end:,}")
        self.logger.info(f"Expansion: -{upstream_bp:,}bp/+{downstream_bp:,}bp")
        self.logger.info(f"Final region: {region_start:,}-{region_end:,}")
        
        return region_start, region_end
    
    def generate_landscape_view(self) -> bool:
        """Render the landscape view according to config."""
        try:
            # Gather parameters
            gene_names = self.config['target']['gene_names']
            value_threshold = self.config['thresholds']['value_threshold']
            
            # Compute genomic window
            region_start, region_end = self.calculate_region()
            
            # Use gene metadata to determine chromosome
            gene_info = self.bulk_api.get_gene_info(gene_names)
            chromosome = gene_info.iloc[0]['Chromosome'] if not gene_info.empty else 'chr1'
            
            self.logger.info("🎨 Generating landscape view...")
            self.logger.info(f"   Genes: {gene_names}")
            self.logger.info(f"   Region: {chromosome}:{region_start:,}-{region_end:,}")
            self.logger.info(f"   Threshold: linkage ≥ {value_threshold}")
            
            # Render Plotly figure
            fig = self.landscape_api.create_landscape_view_plot(
                gene_names=gene_names,
                chromosome=chromosome,
                region_start=region_start,
                region_end=region_end,
                bulk_multi_api=self.bulk_api,
                value_threshold=value_threshold
            )
            
            if not fig or not fig.data:
                self.logger.error("Generated figure is empty")
                return False
            
            self.logger.info(f"📊 Figure generated with {len(fig.data)} tracks, height={fig.layout.height}px")
            
            # Persist to disk
            return self._save_outputs(fig, region_start, region_end)
            
        except Exception as e:
            self.logger.error(f"Failed to generate landscape view: {e}")
            return False
    
    def _save_outputs(self, fig, region_start: int, region_end: int) -> bool:
        """Persist renders to disk."""
        try:
            try:
                import plotly.io as pio
            except ImportError as exc:
                raise RuntimeError("Missing dependency: plotly (install via `pip install plotly`)") from exc

            output_dir = Path(self.config['general']['output_dir'])
            output_name = self.config['general']['output_name']
            formats = self.config['general']['formats']
            
            success_count = 0
            
            # Save HTML
            if 'html' in formats:
                html_config = {
                    'include_plotlyjs': self.config['export'].get('html_include_plotlyjs', True),
                    'config': {'displayModeBar': self.config['export'].get('html_display_modebar', True)},
                    'div_id': self.config['export'].get('html_div_id', 'landscape-view-plot')
                }
                
                html_file = output_dir / f"{output_name}.html"
                pio.write_html(fig, html_file, **html_config)
                self.logger.info(f"💾 Saved HTML: {html_file}")
                success_count += 1
            
            # Save PNG
            if 'png' in formats:
                try:
                    png_file = output_dir / f"{output_name}.png"
                    width = self.config['export']['png_width']
                    height = self.config['export']['png_height']
                    
                    pio.write_image(fig, png_file, width=width, height=height)
                    self.logger.info(f"💾 Saved PNG: {png_file}")
                    success_count += 1
                    
                except Exception as e:
                    self.logger.warning(f"PNG export failed: {e}")
                    self.logger.info("Hint: PNG export requires kaleido (pip install kaleido)")
            
            # Save PDF (experimental)
            if 'pdf' in formats and self.config['export'].get('pdf_enabled', False):
                try:
                    pdf_file = output_dir / f"{output_name}.pdf"
                    width_inch = self.config['export'].get('pdf_width', 14)
                    height_inch = self.config['export'].get('pdf_height', 10)
                    
                    pio.write_image(fig, pdf_file, width=width_inch*72, height=height_inch*72, format='pdf')
                    self.logger.info(f"💾 Saved PDF: {pdf_file}")
                    success_count += 1
                    
                except Exception as e:
                    self.logger.warning(f"PDF export failed: {e}")
            
            return success_count > 0
            
        except Exception as e:
            self.logger.error(f"Failed to save outputs: {e}")
            return False
    


def main():
    """Entry point for CLI execution."""
    parser = argparse.ArgumentParser(
        description='Landscape View CLI - multi-omics landscape renderer',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Run with default config file
  python landscape_view_cli.py

  # Run with a custom config
  python landscape_view_cli.py --config my_config.toml

  # Enable debug logging
  python landscape_view_cli.py --debug
        """
    )
    
    parser.add_argument(
        '--config', '-c',
        type=str,
        default='landscape_view_config.toml',
        help='Path to config file (default: landscape_view_config.toml)'
    )
    
    parser.add_argument(
        '--debug', '-d',
        action='store_true',
        help='Enable verbose logging for debugging'
    )
    
    parser.add_argument(
        '--version', '-v',
        action='version',
        version='Landscape View CLI v1.0.0'
    )
    
    args = parser.parse_args()
    
    # Instantiate CLI
    cli = LandscapeViewCLI()
    
    # Enable debug logging if requested
    if args.debug:
        cli._setup_logging(debug=True)
    
    try:
        # Load config and generate the view
        if not cli.load_config(args.config):
            sys.exit(1)
        
        if not cli.validate_config():
            sys.exit(1)
        
        if not cli.initialize_apis():
            sys.exit(1)
        
        success = cli.generate_landscape_view()
        
        if success:
            print("🎉 Landscape view generated successfully!")
            sys.exit(0)
        else:
            print("❌ Landscape view generation failed!")
            sys.exit(1)
            
    except KeyboardInterrupt:
        print("\n⚠️  Operation interrupted by user")
        sys.exit(1)
    except Exception as e:
        print(f"❌ CLI execution failed: {e}")
        if args.debug:
            import traceback
            traceback.print_exc()
        sys.exit(1)


if __name__ == "__main__":
    main()
