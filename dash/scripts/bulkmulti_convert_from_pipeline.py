#!/usr/bin/env python3
"""
Convert miniENCODE Pipeline outputs to BulkMulti format

This script processes ANANSE outputs from miniENCODE pipeline and converts them
to the BulkMulti format used by Dash Landscape View.

Input:
- Pipeline directory containing ANANSE outputs, ELS peaks, RNA TPM, and BigWig files
- metadata.csv file describing samples
- GTF annotation file

Output:
- BulkMulti.db: SQLite database with 5 tables (genes, exons, peaks, linkages, rna_expression)
- atac_signals/*.bw: BigWig signal tracks
- sample_config.tsv: Sample configuration file
"""

import pandas as pd
import numpy as np
from pathlib import Path
import argparse
import sys
import time
import sqlite3
import shutil
from typing import List, Dict, Tuple, Optional, Set

# Add project root to path for imports
PROJECT_ROOT = next(
    (parent for parent in Path(__file__).resolve().parents if (parent / "backend").exists()),
    None,
)
if PROJECT_ROOT and str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from scripts.logging_utils import configure_script_logger, create_print_proxy

_SCRIPT_LOGGER = configure_script_logger(__name__)
print = create_print_proxy(_SCRIPT_LOGGER)


class BulkMultiConverter:
    """Converter for miniENCODE Pipeline outputs to BulkMulti format"""

    def __init__(
        self,
        pipeline_dir: Path,
        metadata_file: Path,
        gtf_file: Path,
        output_dir: Path,
        species: str,
    ):
        self.pipeline_dir = Path(pipeline_dir)
        self.metadata_file = Path(metadata_file)
        self.gtf_file = Path(gtf_file)
        self.output_dir = Path(output_dir)
        self.species = species

        # Output files
        self.db_file = self.output_dir / "BulkMulti.db"
        self.atac_signals_dir = self.output_dir / "atac_signals"
        self.sample_config_file = self.output_dir / "sample_config.tsv"

        # Data holders
        self.metadata_df = None
        self.samples = []
        self.genes_df = None
        self.exons_df = None
        self.peaks_df = None
        self.linkages_df = None
        self.rna_expression_df = None

    def validate_inputs(self) -> bool:
        """Validate input files and directories"""
        print("🔍 Validating inputs...")

        # Check pipeline directory
        if not self.pipeline_dir.exists():
            print(f"❌ Pipeline directory not found: {self.pipeline_dir}")
            return False

        # Check metadata file
        if not self.metadata_file.exists():
            print(f"❌ Metadata file not found: {self.metadata_file}")
            return False

        # Check GTF file
        if not self.gtf_file.exists():
            print(f"❌ GTF file not found: {self.gtf_file}")
            return False

        print(f"✅ Pipeline directory: {self.pipeline_dir}")
        print(f"✅ Metadata file: {self.metadata_file}")
        print(f"✅ GTF file: {self.gtf_file}")

        return True

    def load_metadata(self) -> bool:
        """Load and validate sample metadata"""
        print("\n📋 Loading sample metadata...")

        try:
            # Read metadata CSV
            self.metadata_df = pd.read_csv(
                self.metadata_file,
                comment='#',
                skipinitialspace=True
            )

            # Required columns
            required_cols = ['sample_id', 'display_name', 'color', 'description']
            missing_cols = [col for col in required_cols if col not in self.metadata_df.columns]
            if missing_cols:
                print(f"❌ Missing required columns: {', '.join(missing_cols)}")
                print(f"   Required: {', '.join(required_cols)}")
                return False

            # Get sample list
            self.samples = self.metadata_df['sample_id'].tolist()
            print(f"✅ Loaded metadata for {len(self.samples)} samples:")
            for sample in self.samples:
                print(f"   - {sample}")

            return True

        except Exception as e:
            print(f"❌ Error loading metadata: {e}")
            return False

    def infer_file_paths(self, sample_id: str) -> Dict[str, Path]:
        """Infer file paths for a sample based on Pipeline directory structure"""
        paths = {
            'ananse_network': self.pipeline_dir / f"data/GRN/{sample_id}/{sample_id}.network.txt",
            'els_bed': self.pipeline_dir / f"data/ELS/{sample_id}/ELS_{sample_id}.bed",
            'tpm_file': self.pipeline_dir / f"data/GRN/input/tpm_{sample_id}.txt",
            'bigwig_file': self.pipeline_dir / f"data/alignment/{sample_id}.bam.bw",
        }

        # Allow override from metadata if columns exist
        row = self.metadata_df[self.metadata_df['sample_id'] == sample_id].iloc[0]
        for key in paths.keys():
            if key in row and pd.notna(row[key]):
                paths[key] = Path(row[key])

        return paths

    def verify_sample_files(self) -> bool:
        """Verify that all required files exist for each sample"""
        print("\n🔍 Verifying sample files...")

        all_valid = True
        for sample_id in self.samples:
            paths = self.infer_file_paths(sample_id)

            print(f"\n  Sample: {sample_id}")
            for file_type, file_path in paths.items():
                path_source = "default"
                row = self.metadata_df[self.metadata_df['sample_id'] == sample_id].iloc[0]
                if file_type in row and pd.notna(row[file_type]):
                    path_source = "metadata"
                if file_path.exists():
                    print(f"    ✅ {file_type} ({path_source}): {file_path.name}")
                else:
                    print(f"    ❌ {file_type} ({path_source}): NOT FOUND - {file_path}")
                    all_valid = False

        if not all_valid:
            print("\n❌ Some required files are missing")
            return False

        print("\n✅ All sample files verified")
        return True

    def parse_gtf(self) -> bool:
        """Parse GTF file to extract genes and exons"""
        print("\n🧬 Parsing GTF annotation...")

        try:
            genes_list = []
            exons_list = []
            gene_id_map = {}  # Map gene_name to auto-increment gene_id
            next_gene_id = 1

            with open(self.gtf_file, 'r') as f:
                for line in f:
                    if line.startswith('#'):
                        continue

                    fields = line.strip().split('\t')
                    if len(fields) < 9:
                        continue

                    chrom, source, feature, start, end, score, strand, frame, attributes = fields

                    # Parse attributes
                    attr_dict = {}
                    for attr in attributes.strip().rstrip(';').split(';'):
                        attr = attr.strip()
                        if ' ' in attr:
                            key, value = attr.split(' ', 1)
                            attr_dict[key] = value.strip('"')

                    gene_name = attr_dict.get('gene_name', attr_dict.get('gene_id', ''))

                    if not gene_name:
                        continue

                    # Process gene feature
                    if feature == 'gene':
                        # Assign gene_id if not seen before
                        if gene_name not in gene_id_map:
                            gene_id_map[gene_name] = next_gene_id
                            next_gene_id += 1

                        gene_id = gene_id_map[gene_name]

                        genes_list.append({
                            'Symbol': gene_name,
                            'Chromosome': chrom,
                            'Start': int(start),
                            'End': int(end),
                            'Strand': strand,
                            'gene_id': gene_id
                        })

                    # Process exon feature
                    elif feature == 'exon':
                        # Assign gene_id if not seen before
                        if gene_name not in gene_id_map:
                            gene_id_map[gene_name] = next_gene_id
                            next_gene_id += 1

                        gene_id = gene_id_map[gene_name]

                        exons_list.append({
                            'Gene_Symbol': gene_name,
                            'Chromosome': chrom,
                            'Start': int(start),
                            'End': int(end),
                            'Strand': strand,
                            'gene_id': gene_id
                        })

            self.genes_df = pd.DataFrame(genes_list)
            self.exons_df = pd.DataFrame(exons_list)

            print(f"✅ Parsed {len(self.genes_df)} genes")
            print(f"✅ Parsed {len(self.exons_df)} exons")

            return True

        except Exception as e:
            print(f"❌ Error parsing GTF: {e}")
            import traceback
            traceback.print_exc()
            return False

    def load_els_peaks(self) -> bool:
        """Load ELS peaks from all samples"""
        print("\n🏔️  Loading ELS peaks...")

        try:
            all_peaks = []
            peak_id_counter = 1

            for sample_id in self.samples:
                paths = self.infer_file_paths(sample_id)
                els_bed = paths['els_bed']

                # Read BED file
                df = pd.read_csv(
                    els_bed,
                    sep='\t',
                    header=None,
                    names=['Chromosome', 'Start', 'End', 'Name', 'Score', 'Strand'],
                    usecols=[0, 1, 2, 3, 4, 5] if els_bed.open('r').readline().count('\t') >= 5 else [0, 1, 2]
                )

                # Standardize columns
                if 'Strand' not in df.columns:
                    df['Strand'] = '*'

                # Add peak_id (this will be used to map ANANSE source)
                df['peak_id'] = [f"peak_{peak_id_counter + i}" for i in range(len(df))]
                peak_id_counter += len(df)

                all_peaks.append(df[['Chromosome', 'Start', 'End', 'Strand', 'peak_id']])
                print(f"  Loaded {len(df)} peaks from {sample_id}")

            # Combine and deduplicate peaks
            self.peaks_df = pd.concat(all_peaks, ignore_index=True)

            # Remove duplicates based on coordinates
            before_dedup = len(self.peaks_df)
            self.peaks_df = self.peaks_df.drop_duplicates(
                subset=['Chromosome', 'Start', 'End'],
                keep='first'
            ).reset_index(drop=True)

            print(f"✅ Loaded {len(self.peaks_df)} unique peaks (deduplicated from {before_dedup})")

            return True

        except Exception as e:
            print(f"❌ Error loading ELS peaks: {e}")
            import traceback
            traceback.print_exc()
            return False

    def load_ananse_networks(self) -> bool:
        """Load ANANSE network files and create linkages table"""
        print("\n🔗 Loading ANANSE networks...")

        try:
            all_linkages = []

            for sample_id in self.samples:
                paths = self.infer_file_paths(sample_id)
                network_file = paths['ananse_network']

                # Read ANANSE network.txt
                df = pd.read_csv(network_file, sep='\t')

                # Expected columns: source, target, binding, factor_expression, target_expression, prob
                required_cols = ['source', 'target', 'prob']
                if not all(col in df.columns for col in required_cols):
                    print(f"❌ Missing required columns in {network_file}")
                    print(f"   Expected: {required_cols}")
                    print(f"   Found: {list(df.columns)}")
                    return False

                print(f"  Loaded {len(df)} connections from {sample_id}")
                all_linkages.append(df[required_cols])

            # Combine all linkages
            combined_linkages = pd.concat(all_linkages, ignore_index=True)

            # Deduplicate by (source, target) pair, keeping highest prob
            combined_linkages = combined_linkages.sort_values('prob', ascending=False)
            combined_linkages = combined_linkages.drop_duplicates(
                subset=['source', 'target'],
                keep='first'
            ).reset_index(drop=True)

            print(f"✅ Loaded {len(combined_linkages)} unique linkages")

            # Map peak source to coordinates and calculate distance
            if not self._map_peaks_and_calculate_distance(combined_linkages):
                return False

            return True

        except Exception as e:
            print(f"❌ Error loading ANANSE networks: {e}")
            import traceback
            traceback.print_exc()
            return False

    def _map_peaks_and_calculate_distance(self, linkages_df: pd.DataFrame) -> bool:
        """Map peak IDs to coordinates and calculate distance to gene TSS"""
        print("\n📏 Mapping peaks and calculating distances...")

        try:
            # Create peak coordinate lookup from peaks_df
            # ANANSE source field contains peak coordinates or IDs
            # We need to parse the source field to extract coordinates

            linkages_list = []

            for idx, row in linkages_df.iterrows():
                source = row['source']
                target = row['target']
                prob = row['prob']

                # Parse source to extract coordinates
                # ANANSE source format is typically: chr:start-end or peak_id
                if ':' in source and '-' in source:
                    # Format: chr:start-end
                    chrom, pos_range = source.split(':', 1)
                    start, end = pos_range.split('-')
                    chrom = chrom.strip()
                    start = int(start)
                    end = int(end)
                else:
                    # Try to find matching peak by ID or skip
                    print(f"⚠️  Warning: Cannot parse source '{source}', skipping")
                    continue

                # Find target gene in genes_df
                gene_rows = self.genes_df[self.genes_df['Symbol'] == target]
                if len(gene_rows) == 0:
                    # Gene not in annotation, skip
                    continue

                gene_row = gene_rows.iloc[0]
                gene_chrom = gene_row['Chromosome']
                gene_start = gene_row['Start']
                gene_end = gene_row['End']
                gene_strand = gene_row['Strand']

                # Skip if chromosomes don't match
                if chrom != gene_chrom:
                    continue

                # Calculate TSS based on strand
                if gene_strand == '+':
                    tss = gene_start
                elif gene_strand == '-':
                    tss = gene_end
                else:
                    tss = gene_start  # Default to start

                # Calculate distance from peak midpoint to TSS
                peak_midpoint = (start + end) // 2
                distance = peak_midpoint - tss

                # Create linkage record
                linkages_list.append({
                    'seqnames': chrom,
                    'start': start,
                    'end': end,
                    'width': end - start,
                    'strand': '*',
                    'gene_name': target,
                    'dist': distance,
                    'value': prob
                })

                if (idx + 1) % 10000 == 0:
                    print(f"  Processed {idx + 1:,} linkages...")

            self.linkages_df = pd.DataFrame(linkages_list)
            print(f"✅ Created {len(self.linkages_df)} linkages with coordinates and distances")

            return True

        except Exception as e:
            print(f"❌ Error mapping peaks and calculating distance: {e}")
            import traceback
            traceback.print_exc()
            return False

    def load_rna_expression(self) -> bool:
        """Load RNA expression data in wide format"""
        print("\n📊 Loading RNA expression data...")

        try:
            # Collect TPM data from all samples
            sample_tpms = {}

            for sample_id in self.samples:
                paths = self.infer_file_paths(sample_id)
                tpm_file = paths['tpm_file']

                # Read TPM file (typically gene_name\tTPM format)
                df = pd.read_csv(tpm_file, sep='\t', header=None, names=['gene_name', 'TPM'])
                sample_tpms[sample_id] = df.set_index('gene_name')['TPM']

                print(f"  Loaded TPM for {len(df)} genes from {sample_id}")

            # Create wide format DataFrame
            rna_wide = pd.DataFrame(sample_tpms)
            rna_wide = rna_wide.reset_index()
            rna_wide = rna_wide.rename(columns={'index': 'name'})

            # Reorder columns to match sample order in metadata
            sample_cols = self.metadata_df['sample_id'].tolist()
            final_cols = sample_cols + ['name']
            self.rna_expression_df = rna_wide[final_cols]

            print(f"✅ Created RNA expression matrix: {len(self.rna_expression_df)} genes × {len(sample_cols)} samples")

            return True

        except Exception as e:
            print(f"❌ Error loading RNA expression: {e}")
            import traceback
            traceback.print_exc()
            return False

    def save_to_sqlite(self) -> bool:
        """Save all tables to SQLite database"""
        print("\n💾 Saving to SQLite database...")

        try:
            # Create output directory
            self.output_dir.mkdir(parents=True, exist_ok=True)

            # Connect to database
            conn = sqlite3.connect(self.db_file)

            # Save tables
            # Note: peaks table doesn't include peak_id in final output
            peaks_output = self.peaks_df[['Chromosome', 'Start', 'End', 'Strand']].copy()

            self.genes_df.to_sql('genes', conn, if_exists='replace', index=False)
            print(f"  ✅ Saved genes table: {len(self.genes_df)} records")

            self.exons_df.to_sql('exons', conn, if_exists='replace', index=False)
            print(f"  ✅ Saved exons table: {len(self.exons_df)} records")

            peaks_output.to_sql('peaks', conn, if_exists='replace', index=False)
            print(f"  ✅ Saved peaks table: {len(peaks_output)} records")

            self.linkages_df.to_sql('linkages', conn, if_exists='replace', index=False)
            print(f"  ✅ Saved linkages table: {len(self.linkages_df)} records")

            self.rna_expression_df.to_sql('rna_expression', conn, if_exists='replace', index=False)
            print(f"  ✅ Saved rna_expression table: {len(self.rna_expression_df)} genes × {len(self.samples)} samples")

            conn.close()
            print(f"✅ Database saved to: {self.db_file}")

            return True

        except Exception as e:
            print(f"❌ Error saving to SQLite: {e}")
            import traceback
            traceback.print_exc()
            return False

    def copy_bigwig_files(self) -> bool:
        """Copy BigWig files to atac_signals directory"""
        print("\n📁 Copying BigWig files...")

        try:
            # Create atac_signals directory
            self.atac_signals_dir.mkdir(parents=True, exist_ok=True)

            for sample_id in self.samples:
                paths = self.infer_file_paths(sample_id)
                bigwig_src = paths['bigwig_file']
                bigwig_dst = self.atac_signals_dir / f"{sample_id}.bw"

                # Copy or create symlink
                if bigwig_src.exists():
                    shutil.copy2(bigwig_src, bigwig_dst)
                    print(f"  ✅ Copied {sample_id}.bw")
                else:
                    print(f"  ⚠️  Warning: BigWig not found for {sample_id}: {bigwig_src}")

            print(f"✅ BigWig files copied to: {self.atac_signals_dir}")

            return True

        except Exception as e:
            print(f"❌ Error copying BigWig files: {e}")
            import traceback
            traceback.print_exc()
            return False

    def generate_sample_config(self) -> bool:
        """Generate sample_config.tsv"""
        print("\n📝 Generating sample config...")

        try:
            # Create sample config from metadata
            config_df = self.metadata_df[['sample_id', 'display_name', 'color', 'description']].copy()
            config_df.columns = ['sample_name', 'display_name', 'color', 'description']

            # Save to TSV
            config_df.to_csv(self.sample_config_file, sep='\t', index=False)

            print(f"✅ Sample config saved to: {self.sample_config_file}")
            print(f"   {len(config_df)} samples configured")

            return True

        except Exception as e:
            print(f"❌ Error generating sample config: {e}")
            import traceback
            traceback.print_exc()
            return False

    def verify_output(self) -> bool:
        """Verify output files and database structure"""
        print("\n✅ Verifying output...")

        try:
            # Check database file
            if not self.db_file.exists():
                print(f"❌ Database file not found: {self.db_file}")
                return False

            # Check database tables
            conn = sqlite3.connect(self.db_file)
            cursor = conn.cursor()

            expected_tables = ['genes', 'exons', 'peaks', 'linkages', 'rna_expression']
            cursor.execute("SELECT name FROM sqlite_master WHERE type='table'")
            actual_tables = [row[0] for row in cursor.fetchall()]

            for table in expected_tables:
                if table not in actual_tables:
                    print(f"❌ Missing table: {table}")
                    conn.close()
                    return False

                cursor.execute(f"SELECT COUNT(*) FROM {table}")
                count = cursor.fetchone()[0]
                print(f"  ✅ {table}: {count:,} records")

            conn.close()

            # Check sample config
            if not self.sample_config_file.exists():
                print(f"❌ Sample config not found: {self.sample_config_file}")
                return False

            print(f"  ✅ sample_config.tsv: {len(self.samples)} samples")

            # Check atac_signals directory
            if not self.atac_signals_dir.exists():
                print(f"❌ atac_signals directory not found: {self.atac_signals_dir}")
                return False

            bigwig_count = len(list(self.atac_signals_dir.glob("*.bw")))
            print(f"  ✅ atac_signals/: {bigwig_count} BigWig files")

            print("\n✅ All output files verified")

            return True

        except Exception as e:
            print(f"❌ Verification error: {e}")
            import traceback
            traceback.print_exc()
            return False

    def run(self, dry_run: bool = False) -> bool:
        """Run the complete conversion pipeline"""
        start_time = time.time()

        print("=" * 70)
        print("🚀 BulkMulti Data Conversion from miniENCODE Pipeline")
        print("=" * 70)
        print(f"Pipeline directory: {self.pipeline_dir}")
        print(f"Metadata file:      {self.metadata_file}")
        print(f"GTF file:           {self.gtf_file}")
        print(f"Output directory:   {self.output_dir}")
        print(f"Species:            {self.species}")
        print("=" * 70)
        print()

        # Step 1: Validate inputs
        if not self.validate_inputs():
            return False

        # Step 2: Load metadata
        if not self.load_metadata():
            return False

        # Step 3: Verify sample files
        if not self.verify_sample_files():
            return False

        # Step 4: Parse GTF
        if not self.parse_gtf():
            return False

        # Step 5: Load ELS peaks
        if not self.load_els_peaks():
            return False

        # Step 6: Load ANANSE networks
        if not self.load_ananse_networks():
            return False

        # Step 7: Load RNA expression
        if not self.load_rna_expression():
            return False

        if dry_run:
            print("\n🔍 Dry run mode - stopping before writing output")
            print(f"   Would generate BulkMulti.db with {len(self.linkages_df):,} linkages")
            print(f"   Output directory: {self.output_dir}")
            return True

        # Step 8: Save to SQLite
        if not self.save_to_sqlite():
            return False

        # Step 9: Copy BigWig files
        if not self.copy_bigwig_files():
            return False

        # Step 10: Generate sample config
        if not self.generate_sample_config():
            return False

        # Step 11: Verify output
        if not self.verify_output():
            return False

        elapsed = time.time() - start_time
        print(f"\n🎉 Conversion completed successfully in {elapsed:.1f}s")
        print(f"📁 Output directory: {self.output_dir}/")

        return True


def main():
    parser = argparse.ArgumentParser(
        description="Convert miniENCODE Pipeline outputs to Dash BulkMulti format",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Basic usage
  python bulkmulti_convert_from_pipeline.py \\
      --species zebrafish \\
      --metadata bulkmulti_metadata.csv \\
      --pipeline-dir /path/to/pipeline \\
      --gtf /path/to/pipeline/data/reference/GRCz11.gtf \\
      --output-dir data/zebrafish/BulkMulti/2025Study

  # Dry run to test without writing
  python bulkmulti_convert_from_pipeline.py \\
      --species zebrafish \\
      --metadata bulkmulti_metadata.csv \\
      --pipeline-dir /path/to/pipeline \\
      --gtf /path/to/pipeline/data/reference/GRCz11.gtf \\
      --output-dir data/zebrafish/BulkMulti/2025Study \\
      --dry-run
"""
    )

    parser.add_argument(
        '--species',
        required=True,
        help='Species identifier (e.g., zebrafish, medaka)'
    )
    parser.add_argument(
        '--metadata',
        required=True,
        help='CSV file with sample metadata (columns: sample_id, display_name, color, description)'
    )
    parser.add_argument(
        '--pipeline-dir',
        required=True,
        help='miniENCODE Pipeline root directory'
    )
    parser.add_argument(
        '--gtf',
        required=True,
        help='Gene annotation GTF file'
    )
    parser.add_argument(
        '--output-dir',
        required=True,
        help='Output directory for BulkMulti data'
    )
    parser.add_argument(
        '--dry-run',
        action='store_true',
        help='Run validation and processing but do not write output files'
    )

    args = parser.parse_args()

    # Create converter
    converter = BulkMultiConverter(
        pipeline_dir=args.pipeline_dir,
        metadata_file=args.metadata,
        gtf_file=args.gtf,
        output_dir=args.output_dir,
        species=args.species,
    )

    # Run conversion
    success = converter.run(dry_run=args.dry_run)

    sys.exit(0 if success else 1)


if __name__ == '__main__':
    main()
