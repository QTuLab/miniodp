#!/usr/bin/env python3
"""
Convert StringTie gene_abundance.tab files to BulkRNA Parquet format

This script processes RNA-seq quantification results from miniENCODE pipeline
(or any StringTie output) and converts them to the standardized Parquet format
used by Dash ExpressionAPI.

Input:
- Directory containing *_gene_abundance.tab files (StringTie output)
- metadata.csv file describing samples

Output:
- expression_data.parquet: Main expression matrix
- gene_index.parquet: Per-gene statistics
- source_index.parquet: Per-source statistics

Notes:
- If multiple rows share the same public Source/Sample pair, ExpressionAPI
  will aggregate them at query time with mean TPM.
- An optional sample_id column is preserved in expression_data.parquet for
  internal traceability.
"""

import pandas as pd
import numpy as np
from pathlib import Path
import argparse
import sys
import time
import re
import sqlite3
from typing import List, Dict, Tuple, Optional

# Add project root to path for imports
PROJECT_ROOT = next(
    (parent for parent in Path(__file__).resolve().parents if (parent / "backend").exists()),
    None,
)
if PROJECT_ROOT and str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from backend.species_adapters import GeneIdAdapter, MedakaAdapter, create_species_adapter
from scripts.logging_utils import configure_script_logger, create_print_proxy

_SCRIPT_LOGGER = configure_script_logger(__name__)
print = create_print_proxy(_SCRIPT_LOGGER)


class BulkRNAConverter:
    """Converter for StringTie output to Dash BulkRNA Parquet format"""

    def __init__(
        self,
        input_dir: Path,
        metadata_file: Path,
        output_dir: Path,
        species: str,
        min_tpm: float = 0.0,
        validate_genes: bool = False,
    ):
        self.input_dir = Path(input_dir)
        self.metadata_file = Path(metadata_file)
        self.output_dir = Path(output_dir)
        self.species = species
        self.min_tpm = min_tpm
        self.validate_genes = validate_genes
        self.adapter = create_species_adapter(species)

        # Output files
        self.expression_file = self.output_dir / "expression_data.parquet"
        self.gene_index_file = self.output_dir / "gene_index.parquet"
        self.source_index_file = self.output_dir / "source_index.parquet"

        # Data holders
        self.metadata_df = None
        self.abundance_files = []
        self.expression_data = None
        self.valid_genes = None
        self.gene_name_to_id = None
        self.gene_name_lower_to_id = None

        # Detect the common ID column used by this species adapter.
        self.gene_id_col = self.adapter.get_common_id_column()

    def _resolve_geneinfo_db(self) -> Path:
        """Resolve the geneinfo.db path for the selected species."""
        candidates = [
            self.output_dir.parent / "geneinfo.db",
            PROJECT_ROOT / "data" / self.species / "geneinfo.db" if PROJECT_ROOT else None,
        ]
        for candidate in candidates:
            if candidate and candidate.exists():
                return candidate
        return candidates[0]

    def validate_inputs(self) -> bool:
        """Validate input directory and metadata file"""
        print("🔍 Validating inputs...")

        # Check input directory
        if not self.input_dir.exists():
            print(f"❌ Input directory not found: {self.input_dir}")
            return False

        # Find gene_abundance.tab files
        self.abundance_files = sorted(self.input_dir.glob("*_gene_abundance.tab"))
        if not self.abundance_files:
            print(f"❌ No *_gene_abundance.tab files found in {self.input_dir}")
            print("   Expected StringTie output files")
            return False

        print(f"✅ Found {len(self.abundance_files)} gene_abundance.tab files")

        # Check metadata file
        if not self.metadata_file.exists():
            print(f"❌ Metadata file not found: {self.metadata_file}")
            return False

        print(f"✅ Found metadata file: {self.metadata_file.name}")

        return True

    def load_metadata(self) -> bool:
        """Load and validate metadata CSV"""
        print("📋 Loading metadata...")

        try:
            self.metadata_df = pd.read_csv(self.metadata_file)

            # Check required columns
            required_cols = ['sample_id', 'source', 'sample_name']
            missing_cols = [col for col in required_cols if col not in self.metadata_df.columns]
            if missing_cols:
                print(f"❌ Missing required columns in metadata: {missing_cols}")
                print(f"   Found columns: {list(self.metadata_df.columns)}")
                return False

            # Remove any whitespace from strings
            for col in ['sample_id', 'source', 'sample_name']:
                self.metadata_df[col] = self.metadata_df[col].astype(str).str.strip()

            print(f"✅ Loaded {len(self.metadata_df)} sample entries")
            print(f"   - Unique sources: {self.metadata_df['source'].nunique()}")
            print(f"   - Sample IDs: {', '.join(self.metadata_df['sample_id'].head(5).tolist())}")
            if len(self.metadata_df) > 5:
                print(f"     ... and {len(self.metadata_df) - 5} more")

            return True

        except Exception as e:
            print(f"❌ Error loading metadata: {e}")
            return False

    def _extract_sample_id(self, filename: str) -> Optional[str]:
        """Extract sample ID from gene_abundance.tab filename"""
        # Pattern: {sample_id}_gene_abundance.tab
        match = re.match(r"(.+?)_gene_abundance\.tab", filename)
        if match:
            return match.group(1)
        return None

    def _extract_common_id(self, gene_id: str) -> Optional[str]:
        """Extract the species-common gene ID from a StringTie Gene ID field."""
        if pd.isna(gene_id):
            return None

        gene_id = str(gene_id).strip()
        if not gene_id:
            return None

        if isinstance(self.adapter, MedakaAdapter):
            if gene_id.startswith("gene:"):
                gene_id = gene_id.split(":", 1)[1]
            if gene_id.startswith("IGDB"):
                return gene_id.split('.')[0]
            return None

        if isinstance(self.adapter, GeneIdAdapter):
            if gene_id.isdigit():
                return gene_id

            match = re.fullmatch(r"LOC(\d+)", gene_id)
            if match:
                return match.group(1)

            if self.gene_name_to_id and gene_id in self.gene_name_to_id:
                return self.gene_name_to_id[gene_id]

            lowered = gene_id.lower()
            if self.gene_name_lower_to_id and lowered in self.gene_name_lower_to_id:
                return self.gene_name_lower_to_id[lowered]

            match = re.search(r'GeneID:(\d+)', gene_id)
            if match:
                return match.group(1)

            return None

        match = re.search(r'(ENS[A-Z]*[0-9]+)', gene_id)
        if match:
            return match.group(1)

        return None

    def load_gene_database(self) -> bool:
        """Load valid gene IDs from geneinfo.db for validation"""
        geneinfo_db = self._resolve_geneinfo_db()
        need_lookup_maps = isinstance(self.adapter, GeneIdAdapter)

        if not geneinfo_db.exists():
            if need_lookup_maps:
                print(f"❌ Required geneinfo.db not found for {self.species}: {geneinfo_db}")
                return False
            if not self.validate_genes:
                return True
            print(f"⚠️  Gene validation requested but geneinfo.db not found: {geneinfo_db}")
            print("   Proceeding without validation...")
            self.validate_genes = False
            return True

        if not self.validate_genes and not need_lookup_maps:
            return True

        print(f"🧬 Loading gene database for validation...")

        try:
            with sqlite3.connect(geneinfo_db) as conn:
                if need_lookup_maps:
                    gene_df = pd.read_sql_query(
                        "SELECT ENS, GeneName FROM ENS_INFO",
                        conn,
                    )
                    gene_df["ENS"] = gene_df["ENS"].astype(str).str.strip()
                    gene_df["GeneName"] = gene_df["GeneName"].fillna("").astype(str).str.strip()

                    exact_map: Dict[str, str] = {}
                    exact_counts: Dict[str, set[str]] = {}
                    lower_counts: Dict[str, set[str]] = {}

                    for _, row in gene_df.iterrows():
                        gene_id = row["ENS"]
                        gene_name = row["GeneName"]
                        if not gene_id or not gene_name:
                            continue
                        exact_counts.setdefault(gene_name, set()).add(gene_id)
                        lower_counts.setdefault(gene_name.lower(), set()).add(gene_id)

                    self.gene_name_to_id = {
                        gene_name: next(iter(gene_ids))
                        for gene_name, gene_ids in exact_counts.items()
                        if len(gene_ids) == 1
                    }
                    self.gene_name_lower_to_id = {
                        gene_name: next(iter(gene_ids))
                        for gene_name, gene_ids in lower_counts.items()
                        if len(gene_ids) == 1
                    }
                    print(
                        "✅ Loaded killifish GeneName lookup maps: "
                        f"exact={len(self.gene_name_to_id):,}, lower={len(self.gene_name_lower_to_id):,}"
                    )

                if self.validate_genes:
                    query = f"SELECT DISTINCT {self.gene_id_col} FROM ENS_INFO"
                    df = pd.read_sql_query(query, conn)
                    self.valid_genes = set(df[self.gene_id_col].dropna().astype(str))
                    print(f"✅ Loaded {len(self.valid_genes)} valid gene IDs from geneinfo.db")

                return True

        except Exception as e:
            print(f"⚠️  Error loading gene database: {e}")
            if need_lookup_maps:
                return False
            print("   Proceeding without validation...")
            self.validate_genes = False
            return True

    def read_stringtie_outputs(self) -> bool:
        """Read all StringTie gene_abundance.tab files"""
        print(f"📊 Reading {len(self.abundance_files)} StringTie output files...")

        all_data = []
        skipped_samples = []

        for i, file_path in enumerate(self.abundance_files, 1):
            try:
                # Extract sample_id from filename
                sample_id = self._extract_sample_id(file_path.name)
                if not sample_id:
                    print(f"⚠️  Could not extract sample_id from: {file_path.name}")
                    skipped_samples.append(file_path.name)
                    continue

                # Look up metadata
                metadata_row = self.metadata_df[self.metadata_df['sample_id'] == sample_id]
                if metadata_row.empty:
                    print(f"⚠️  No metadata found for sample_id: {sample_id}")
                    skipped_samples.append(sample_id)
                    continue

                source = metadata_row.iloc[0]['source']
                sample_name = metadata_row.iloc[0]['sample_name']
                metadata_sample_id = metadata_row.iloc[0]['sample_id']

                if metadata_sample_id != sample_id:
                    print(
                        "⚠️  Metadata sample_id does not match filename: "
                        f"{metadata_sample_id} vs {sample_id}"
                    )

                # Read StringTie output
                # Columns: Gene ID, Gene Name, Reference, Strand, Start, End, Coverage, FPKM, TPM
                df = pd.read_csv(file_path, sep='\t')

                # Extract species-common IDs used by Dash storage.
                df['gene_id_clean'] = df['Gene ID'].apply(self._extract_common_id)

                # Filter out rows without a valid species-common gene ID
                df_valid = df[df['gene_id_clean'].notna()].copy()

                if len(df_valid) == 0:
                    print(f"⚠️  No valid gene IDs found in: {file_path.name}")
                    skipped_samples.append(file_path.name)
                    continue

                # Build expression records
                df_valid['Source'] = source
                df_valid['Sample'] = sample_name
                df_valid['sample_id'] = metadata_sample_id
                df_valid['Sample_Source'] = source + '::' + sample_name

                # Select and rename columns
                df_expr = df_valid[
                    ['gene_id_clean', 'Source', 'Sample', 'sample_id', 'TPM', 'Sample_Source']
                ].copy()
                df_expr.columns = [
                    self.gene_id_col,
                    'Source',
                    'Sample',
                    'sample_id',
                    'TPM',
                    'Sample_Source',
                ]

                all_data.append(df_expr)

                if i % 10 == 0:
                    print(f"   Processed {i}/{len(self.abundance_files)} files...")

            except Exception as e:
                print(f"⚠️  Error reading {file_path.name}: {e}")
                skipped_samples.append(file_path.name)
                continue

        if not all_data:
            print("❌ No valid data loaded from StringTie files")
            return False

        # Concatenate all data
        self.expression_data = pd.concat(all_data, ignore_index=True)

        key_cols = [self.gene_id_col, 'Source', 'Sample', 'sample_id']
        duplicate_count = int(self.expression_data.duplicated(key_cols).sum())
        if duplicate_count:
            print(
                "⚠️  Found duplicate gene/source/sample/sample_id records: "
                f"{duplicate_count:,}; aggregating TPM by mean"
            )
            self.expression_data = (
                self.expression_data
                .groupby(key_cols, as_index=False, sort=False)
                .agg({'TPM': 'mean'})
            )
            self.expression_data['Sample_Source'] = (
                self.expression_data['Source'] + '::' + self.expression_data['Sample']
            )
            self.expression_data = self.expression_data[
                [self.gene_id_col, 'Source', 'Sample', 'sample_id', 'TPM', 'Sample_Source']
            ]

        print(f"✅ Loaded {len(self.expression_data):,} expression records")
        print(f"   - Unique genes: {self.expression_data[self.gene_id_col].nunique():,}")
        print(f"   - Unique sources: {self.expression_data['Source'].nunique()}")
        print(f"   - Unique samples: {self.expression_data['Sample'].nunique()}")

        if skipped_samples:
            print(f"⚠️  Skipped {len(skipped_samples)} samples due to errors")

        return True

    def validate_gene_ids(self) -> bool:
        """Validate gene IDs against geneinfo.db"""
        if not self.validate_genes or self.valid_genes is None:
            return True

        print("🔍 Validating gene IDs against geneinfo.db...")

        current_genes = set(self.expression_data[self.gene_id_col].unique())
        invalid_genes = current_genes - self.valid_genes

        if invalid_genes:
            print(f"⚠️  Found {len(invalid_genes)} gene IDs not in geneinfo.db:")
            print(f"   Examples: {list(invalid_genes)[:10]}")
            print(f"   These genes will be kept in the output but may not be searchable")
        else:
            print(f"✅ All {len(current_genes):,} genes validated against geneinfo.db")

        return True

    def filter_low_expression(self) -> bool:
        """Filter out low-expression records"""
        if self.min_tpm <= 0:
            return True

        print(f"🔧 Filtering records with TPM < {self.min_tpm}...")

        original_count = len(self.expression_data)
        self.expression_data = self.expression_data[self.expression_data['TPM'] >= self.min_tpm]
        filtered_count = original_count - len(self.expression_data)

        print(f"   Removed {filtered_count:,} low-expression records ({filtered_count/original_count*100:.1f}%)")
        print(f"   Remaining: {len(self.expression_data):,} records")

        return True

    def build_gene_index(self) -> pd.DataFrame:
        """Build gene index for fast lookups"""
        print("🔍 Building gene index...")

        gene_index = self.expression_data.groupby(self.gene_id_col).agg({
            'Source': 'nunique',
            'Sample': 'nunique',
            'TPM': ['mean', 'max', 'count']
        }).round(3)

        # Flatten column names
        gene_index.columns = ['sources', 'samples', 'mean_tpm', 'max_tpm', 'records']

        print(f"✅ Gene index: {len(gene_index):,} genes")

        return gene_index

    def build_source_index(self) -> pd.DataFrame:
        """Build source index for fast lookups"""
        print("🔍 Building source index...")

        source_index = self.expression_data.groupby('Source').agg({
            self.gene_id_col: 'nunique',
            'Sample': 'nunique',
            'TPM': ['mean', 'count']
        }).round(3)

        # Flatten column names
        source_index.columns = ['genes', 'samples', 'mean_tpm', 'records']

        print(f"✅ Source index: {len(source_index)} sources")

        return source_index

    def save_to_parquet(self, gene_index: pd.DataFrame, source_index: pd.DataFrame) -> bool:
        """Save all data to Parquet format"""
        print(f"💾 Saving to Parquet format...")

        # Create output directory
        self.output_dir.mkdir(parents=True, exist_ok=True)

        try:
            # Save expression data
            print(f"   Writing {self.expression_file.name}...")
            self.expression_data.to_parquet(
                self.expression_file,
                compression='snappy',
                index=False
            )

            # Save gene index
            print(f"   Writing {self.gene_index_file.name}...")
            gene_index.to_parquet(self.gene_index_file, compression='snappy')

            # Save source index
            print(f"   Writing {self.source_index_file.name}...")
            source_index.to_parquet(self.source_index_file, compression='snappy')

            # Report file sizes
            expr_size = self.expression_file.stat().st_size / (1024 * 1024)
            gene_size = self.gene_index_file.stat().st_size / (1024 * 1024)
            source_size = self.source_index_file.stat().st_size / (1024 * 1024)
            total_size = expr_size + gene_size + source_size

            print(f"✅ Files saved to {self.output_dir}/")
            print(f"   - expression_data.parquet: {expr_size:.1f} MB")
            print(f"   - gene_index.parquet: {gene_size:.1f} MB")
            print(f"   - source_index.parquet: {source_size:.1f} MB")
            print(f"   Total: {total_size:.1f} MB")

            return True

        except Exception as e:
            print(f"❌ Error saving Parquet files: {e}")
            return False

    def verify_output(self) -> bool:
        """Verify output files"""
        print("🧪 Verifying output...")

        try:
            # Check files exist
            for file_path in [self.expression_file, self.gene_index_file, self.source_index_file]:
                if not file_path.exists():
                    print(f"❌ Missing output file: {file_path}")
                    return False

            # Test loading
            df_expr = pd.read_parquet(self.expression_file)
            df_gene = pd.read_parquet(self.gene_index_file)
            df_source = pd.read_parquet(self.source_index_file)

            print(f"✅ All files verified:")
            print(f"   - Expression data: {len(df_expr):,} records")
            print(f"   - Gene index: {len(df_gene):,} genes")
            print(f"   - Source index: {len(df_source)} sources")

            # Quick stats
            print(f"\n📊 Final statistics:")
            print(f"   - Gene ID column: {self.gene_id_col}")
            print(f"   - TPM range: {df_expr['TPM'].min():.3f} - {df_expr['TPM'].max():.3f}")
            print(f"   - Mean TPM: {df_expr['TPM'].mean():.3f}")
            print(f"   - Median TPM: {df_expr['TPM'].median():.3f}")

            return True

        except Exception as e:
            print(f"❌ Verification error: {e}")
            return False

    def run(self, dry_run: bool = False) -> bool:
        """Run the complete conversion pipeline"""
        start_time = time.time()

        print("=" * 70)
        print("🚀 Bulk RNA-seq Data Conversion")
        print("=" * 70)
        print(f"Input directory: {self.input_dir}")
        print(f"Metadata file:   {self.metadata_file}")
        print(f"Output directory: {self.output_dir}")
        print(f"Species:         {self.species}")
        print(f"Gene ID column:  {self.gene_id_col}")
        if self.min_tpm > 0:
            print(f"Min TPM filter:  {self.min_tpm}")
        print("=" * 70)
        print()

        # Step 1: Validate inputs
        if not self.validate_inputs():
            return False

        # Step 2: Load metadata
        if not self.load_metadata():
            return False

        # Step 3: Load gene database for validation
        if not self.load_gene_database():
            return False

        # Step 4: Read StringTie outputs
        if not self.read_stringtie_outputs():
            return False

        # Step 5: Validate gene IDs
        if not self.validate_gene_ids():
            return False

        # Step 6: Filter low expression
        if not self.filter_low_expression():
            return False

        if dry_run:
            print("\n🔍 Dry run mode - stopping before writing output")
            print(f"   Would generate {len(self.expression_data):,} expression records")
            print(f"   Output directory: {self.output_dir}")
            return True

        # Step 7: Build indexes
        gene_index = self.build_gene_index()
        source_index = self.build_source_index()

        # Step 8: Save to Parquet
        if not self.save_to_parquet(gene_index, source_index):
            return False

        # Step 9: Verify output
        if not self.verify_output():
            return False

        elapsed = time.time() - start_time
        print(f"\n🎉 Conversion completed successfully in {elapsed:.1f}s")
        print(f"📁 Output files: {self.output_dir}/")

        return True


def main():
    parser = argparse.ArgumentParser(
        description="Convert StringTie gene_abundance.tab files to Dash BulkRNA Parquet format",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Basic usage
  python bulkrna_convert_from_stringtie.py \\
      --input-dir /path/to/pipeline/quant \\
      --metadata samples_metadata.csv \\
      --output-dir data/zebrafish/bulkRNA \\
      --species zebrafish

  # With gene validation and TPM filtering
  python bulkrna_convert_from_stringtie.py \\
      --input-dir /path/to/pipeline/quant \\
      --metadata samples_metadata.csv \\
      --output-dir data/zebrafish/bulkRNA \\
      --species zebrafish \\
      --validate-genes \\
      --min-tpm 0.1

  # Dry run to test without writing
  python bulkrna_convert_from_stringtie.py \\
      --input-dir /path/to/pipeline/quant \\
      --metadata samples_metadata.csv \\
      --output-dir data/zebrafish/bulkRNA \\
      --species zebrafish \\
      --dry-run
"""
    )

    parser.add_argument(
        '--input-dir',
        required=True,
        help='Directory containing *_gene_abundance.tab files (StringTie output)'
    )
    parser.add_argument(
        '--metadata',
        required=True,
        help='CSV file with sample metadata (columns: sample_id, source, sample_name)'
    )
    parser.add_argument(
        '--output-dir',
        required=True,
        help='Output directory for bulkRNA Parquet files'
    )
    parser.add_argument(
        '--species',
        required=True,
        help='Species identifier (e.g., zebrafish, medaka)'
    )
    parser.add_argument(
        '--min-tpm',
        type=float,
        default=0.0,
        help='Minimum TPM threshold for filtering (default: 0.0, no filtering)'
    )
    parser.add_argument(
        '--validate-genes',
        action='store_true',
        help='Validate gene IDs against geneinfo.db'
    )
    parser.add_argument(
        '--dry-run',
        action='store_true',
        help='Run validation and processing but do not write output files'
    )

    args = parser.parse_args()

    # Create converter
    converter = BulkRNAConverter(
        input_dir=args.input_dir,
        metadata_file=args.metadata,
        output_dir=args.output_dir,
        species=args.species,
        min_tpm=args.min_tpm,
        validate_genes=args.validate_genes,
    )

    # Run conversion
    success = converter.run(dry_run=args.dry_run)

    return 0 if success else 1


if __name__ == "__main__":
    sys.exit(main())
