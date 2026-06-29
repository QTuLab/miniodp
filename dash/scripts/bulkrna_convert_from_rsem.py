#!/usr/bin/env python3
"""
Convert RSEM genes.results files to Dash BulkRNA Parquet format.

Input:
- Directory containing *.genes.results files
- metadata.csv file describing samples

Output:
- expression_data.parquet: Main expression matrix
- gene_index.parquet: Per-gene statistics
- source_index.parquet: Per-source statistics
"""

from __future__ import annotations

import argparse
import re
import sqlite3
import sys
import time
from pathlib import Path
from typing import Optional

import pandas as pd

# Add project root to path for imports
PROJECT_ROOT = next(
    (parent for parent in Path(__file__).resolve().parents if (parent / "backend").exists()),
    None,
)
if PROJECT_ROOT and str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from backend.species_adapters import MedakaAdapter, create_species_adapter
from scripts.logging_utils import configure_script_logger, create_print_proxy

_SCRIPT_LOGGER = configure_script_logger(__name__)
print = create_print_proxy(_SCRIPT_LOGGER)


class BulkRNARsemConverter:
    """Converter for RSEM gene-level output to Dash BulkRNA Parquet format."""

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

        self.expression_file = self.output_dir / "expression_data.parquet"
        self.gene_index_file = self.output_dir / "gene_index.parquet"
        self.source_index_file = self.output_dir / "source_index.parquet"

        self.metadata_df: Optional[pd.DataFrame] = None
        self.result_files: list[Path] = []
        self.expression_data: Optional[pd.DataFrame] = None
        self.valid_genes: Optional[set[str]] = None

        self.gene_id_col = self.adapter.get_common_id_column()

    def validate_inputs(self) -> bool:
        """Validate input directory and metadata file."""
        print("🔍 Validating inputs...")

        if not self.input_dir.exists():
            print(f"❌ Input directory not found: {self.input_dir}")
            return False

        self.result_files = sorted(self.input_dir.glob("*.genes.results"))
        if not self.result_files:
            print(f"❌ No *.genes.results files found in {self.input_dir}")
            return False

        print(f"✅ Found {len(self.result_files)} genes.results files")

        if not self.metadata_file.exists():
            print(f"❌ Metadata file not found: {self.metadata_file}")
            return False

        print(f"✅ Found metadata file: {self.metadata_file.name}")
        return True

    def load_metadata(self) -> bool:
        """Load metadata CSV while skipping template comment lines."""
        print("📋 Loading metadata...")

        try:
            self.metadata_df = pd.read_csv(self.metadata_file, comment="#")

            required_cols = ["sample_id", "source", "sample_name"]
            missing_cols = [col for col in required_cols if col not in self.metadata_df.columns]
            if missing_cols:
                print(f"❌ Missing required columns in metadata: {missing_cols}")
                print(f"   Found columns: {list(self.metadata_df.columns)}")
                return False

            self.metadata_df = self.metadata_df.dropna(subset=["sample_id", "source", "sample_name"]).copy()
            for col in ["sample_id", "source", "sample_name"]:
                self.metadata_df[col] = self.metadata_df[col].astype(str).str.strip()

            self.metadata_df = self.metadata_df[self.metadata_df["sample_id"] != ""].copy()

            print(f"✅ Loaded {len(self.metadata_df)} sample entries")
            print(f"   - Unique sources: {self.metadata_df['source'].nunique()}")
            print(f"   - Sample IDs: {', '.join(self.metadata_df['sample_id'].head(5).tolist())}")
            if len(self.metadata_df) > 5:
                print(f"     ... and {len(self.metadata_df) - 5} more")

            return True

        except Exception as exc:
            print(f"❌ Error loading metadata: {exc}")
            return False

    def _extract_sample_id(self, filename: str) -> Optional[str]:
        """Extract sample ID from genes.results filename."""
        match = re.match(r"(.+?)\.genes\.results$", filename)
        if match:
            return match.group(1)
        return None

    def _extract_common_id(self, gene_id: str) -> Optional[str]:
        """Normalize a gene ID to the species-common gene identifier."""
        if pd.isna(gene_id):
            return None

        gene_id = str(gene_id).strip()
        if not gene_id:
            return None

        if isinstance(self.adapter, MedakaAdapter):
            if gene_id.startswith("IGDB"):
                return gene_id.split(".", 1)[0]
            return None

        detected = self.adapter.detect_id_format(gene_id)
        if detected in {"ens_id", "igdb_id", "primary_id", "secondary_id"}:
            return gene_id

        match = re.search(r"GeneID:(\d+)", gene_id)
        if match:
            return match.group(1)

        return None

    def load_gene_database(self) -> bool:
        """Load valid gene IDs from geneinfo.db for validation."""
        if not self.validate_genes:
            return True

        geneinfo_db = self.output_dir.parent / "geneinfo.db"
        if not geneinfo_db.exists():
            print(f"⚠️  Gene validation requested but geneinfo.db not found: {geneinfo_db}")
            print("   Proceeding without validation...")
            self.validate_genes = False
            return True

        print("🧬 Loading gene database for validation...")

        try:
            with sqlite3.connect(geneinfo_db) as conn:
                query = f"SELECT DISTINCT {self.gene_id_col} FROM ENS_INFO"
                df = pd.read_sql_query(query, conn)
                self.valid_genes = set(df[self.gene_id_col].dropna())
            print(f"✅ Loaded {len(self.valid_genes)} valid gene IDs from geneinfo.db")
            return True
        except Exception as exc:
            print(f"⚠️  Error loading gene database: {exc}")
            print("   Proceeding without validation...")
            self.validate_genes = False
            return True

    def read_rsem_outputs(self) -> bool:
        """Read all RSEM genes.results files."""
        print(f"📊 Reading {len(self.result_files)} RSEM output files...")

        all_data: list[pd.DataFrame] = []
        skipped_samples: list[str] = []

        for i, file_path in enumerate(self.result_files, 1):
            try:
                sample_id = self._extract_sample_id(file_path.name)
                if not sample_id:
                    print(f"⚠️  Could not extract sample_id from: {file_path.name}")
                    skipped_samples.append(file_path.name)
                    continue

                metadata_row = self.metadata_df[self.metadata_df["sample_id"] == sample_id]
                if metadata_row.empty:
                    print(f"⚠️  No metadata found for sample_id: {sample_id}")
                    skipped_samples.append(sample_id)
                    continue

                source = metadata_row.iloc[0]["source"]
                sample_name = metadata_row.iloc[0]["sample_name"]

                df = pd.read_csv(file_path, sep="\t")
                required_cols = {"gene_id", "TPM"}
                if not required_cols.issubset(df.columns):
                    print(f"⚠️  Missing required columns in {file_path.name}: {sorted(required_cols - set(df.columns))}")
                    skipped_samples.append(file_path.name)
                    continue

                df["gene_id_clean"] = df["gene_id"].apply(self._extract_common_id)
                df_valid = df[df["gene_id_clean"].notna()].copy()

                if df_valid.empty:
                    print(f"⚠️  No valid gene IDs found in: {file_path.name}")
                    skipped_samples.append(file_path.name)
                    continue

                if df_valid["gene_id_clean"].duplicated().any():
                    duplicate_count = int(df_valid["gene_id_clean"].duplicated().sum())
                    print(f"⚠️  Found {duplicate_count} duplicated gene IDs in {file_path.name}; summing TPM by gene")
                    df_valid = df_valid.groupby("gene_id_clean", as_index=False)["TPM"].sum()
                else:
                    df_valid = df_valid[["gene_id_clean", "TPM"]].copy()

                df_valid["Source"] = source
                df_valid["Sample"] = sample_name
                df_valid["Sample_Source"] = source + "::" + sample_name

                df_expr = df_valid[["gene_id_clean", "Source", "Sample", "TPM", "Sample_Source"]].copy()
                df_expr.columns = [self.gene_id_col, "Source", "Sample", "TPM", "Sample_Source"]
                all_data.append(df_expr)

                if i % 10 == 0:
                    print(f"   Processed {i}/{len(self.result_files)} files...")

            except Exception as exc:
                print(f"⚠️  Error reading {file_path.name}: {exc}")
                skipped_samples.append(file_path.name)

        if not all_data:
            print("❌ No valid data loaded from RSEM files")
            return False

        self.expression_data = pd.concat(all_data, ignore_index=True)
        print(f"✅ Loaded {len(self.expression_data):,} expression records")
        print(f"   - Unique genes: {self.expression_data[self.gene_id_col].nunique():,}")
        print(f"   - Unique sources: {self.expression_data['Source'].nunique()}")
        print(f"   - Unique samples: {self.expression_data['Sample'].nunique()}")

        if skipped_samples:
            print(f"⚠️  Skipped {len(skipped_samples)} samples due to errors")

        return True

    def validate_gene_ids(self) -> bool:
        """Validate gene IDs against geneinfo.db."""
        if not self.validate_genes or self.valid_genes is None:
            return True

        print("🔍 Validating gene IDs against geneinfo.db...")
        current_genes = set(self.expression_data[self.gene_id_col].unique())
        invalid_genes = current_genes - self.valid_genes

        if invalid_genes:
            print(f"⚠️  Found {len(invalid_genes)} gene IDs not in geneinfo.db:")
            print(f"   Examples: {list(sorted(invalid_genes))[:10]}")
            print("   These genes will be kept in the output but may not be searchable")
        else:
            print(f"✅ All {len(current_genes):,} genes validated against geneinfo.db")

        return True

    def filter_low_expression(self) -> bool:
        """Filter out low-expression records."""
        if self.min_tpm <= 0:
            return True

        print(f"🔧 Filtering records with TPM < {self.min_tpm}...")
        original_count = len(self.expression_data)
        self.expression_data = self.expression_data[self.expression_data["TPM"] >= self.min_tpm]
        filtered_count = original_count - len(self.expression_data)

        print(f"   Removed {filtered_count:,} low-expression records ({filtered_count / original_count * 100:.1f}%)")
        print(f"   Remaining: {len(self.expression_data):,} records")
        return True

    def build_gene_index(self) -> pd.DataFrame:
        """Build gene index for fast lookups."""
        print("🔍 Building gene index...")
        gene_index = self.expression_data.groupby(self.gene_id_col).agg(
            {
                "Source": "nunique",
                "Sample": "nunique",
                "TPM": ["mean", "max", "count"],
            }
        ).round(3)
        gene_index.columns = ["sources", "samples", "mean_tpm", "max_tpm", "records"]
        print(f"✅ Gene index: {len(gene_index):,} genes")
        return gene_index

    def build_source_index(self) -> pd.DataFrame:
        """Build source index for fast lookups."""
        print("🔍 Building source index...")
        source_index = self.expression_data.groupby("Source").agg(
            {
                self.gene_id_col: "nunique",
                "Sample": "nunique",
                "TPM": ["mean", "count"],
            }
        ).round(3)
        source_index.columns = ["genes", "samples", "mean_tpm", "records"]
        print(f"✅ Source index: {len(source_index)} sources")
        return source_index

    def save_to_parquet(self, gene_index: pd.DataFrame, source_index: pd.DataFrame) -> bool:
        """Save all outputs to Parquet."""
        print("💾 Saving to Parquet format...")
        self.output_dir.mkdir(parents=True, exist_ok=True)

        try:
            print(f"   Writing {self.expression_file.name}...")
            self.expression_data.to_parquet(self.expression_file, compression="snappy", index=False)

            print(f"   Writing {self.gene_index_file.name}...")
            gene_index.to_parquet(self.gene_index_file, compression="snappy")

            print(f"   Writing {self.source_index_file.name}...")
            source_index.to_parquet(self.source_index_file, compression="snappy")

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
        except Exception as exc:
            print(f"❌ Error saving Parquet files: {exc}")
            return False

    def verify_output(self) -> bool:
        """Verify output files."""
        print("🧪 Verifying output...")

        try:
            for file_path in [self.expression_file, self.gene_index_file, self.source_index_file]:
                if not file_path.exists():
                    print(f"❌ Missing output file: {file_path}")
                    return False

            df_expr = pd.read_parquet(self.expression_file)
            df_gene = pd.read_parquet(self.gene_index_file)
            df_source = pd.read_parquet(self.source_index_file)

            print("✅ All files verified:")
            print(f"   - Expression data: {len(df_expr):,} records")
            print(f"   - Gene index: {len(df_gene):,} genes")
            print(f"   - Source index: {len(df_source)} sources")
            print("\n📊 Final statistics:")
            print(f"   - Gene ID column: {self.gene_id_col}")
            print(f"   - TPM range: {df_expr['TPM'].min():.3f} - {df_expr['TPM'].max():.3f}")
            print(f"   - Mean TPM: {df_expr['TPM'].mean():.3f}")
            print(f"   - Median TPM: {df_expr['TPM'].median():.3f}")
            return True
        except Exception as exc:
            print(f"❌ Verification error: {exc}")
            return False

    def run(self, dry_run: bool = False) -> bool:
        """Run the complete conversion pipeline."""
        start_time = time.time()

        print("=" * 70)
        print("🚀 Bulk RNA-seq Data Conversion from RSEM")
        print("=" * 70)
        print(f"Input directory:  {self.input_dir}")
        print(f"Metadata file:    {self.metadata_file}")
        print(f"Output directory: {self.output_dir}")
        print(f"Species:          {self.species}")
        print(f"Gene ID column:   {self.gene_id_col}")
        if self.min_tpm > 0:
            print(f"Min TPM filter:   {self.min_tpm}")
        print("=" * 70)
        print()

        if not self.validate_inputs():
            return False
        if not self.load_metadata():
            return False
        if not self.load_gene_database():
            return False
        if not self.read_rsem_outputs():
            return False
        if not self.validate_gene_ids():
            return False
        if not self.filter_low_expression():
            return False

        if dry_run:
            print("\n🔍 Dry run mode - stopping before writing output")
            print(f"   Would generate {len(self.expression_data):,} expression records")
            print(f"   Output directory: {self.output_dir}")
            return True

        gene_index = self.build_gene_index()
        source_index = self.build_source_index()
        if not self.save_to_parquet(gene_index, source_index):
            return False
        if not self.verify_output():
            return False

        elapsed = time.time() - start_time
        print(f"\n🎉 Conversion completed successfully in {elapsed:.1f}s")
        print(f"📁 Output files: {self.output_dir}/")
        return True


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Convert RSEM genes.results files to Dash BulkRNA Parquet format",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python bulkrna_convert_from_rsem.py \\
      --input-dir /path/to/rsem_results \\
      --metadata metadata.csv \\
      --output-dir data/german_cockroach/bulkRNA \\
      --species german_cockroach

  python bulkrna_convert_from_rsem.py \\
      --input-dir /path/to/rsem_results \\
      --metadata metadata.csv \\
      --output-dir data/german_cockroach/bulkRNA \\
      --species german_cockroach \\
      --validate-genes \\
      --min-tpm 0.1
        """,
    )
    parser.add_argument("--input-dir", required=True, help="Directory containing *.genes.results files")
    parser.add_argument("--metadata", required=True, help="Sample metadata CSV file")
    parser.add_argument("--output-dir", required=True, help="Output directory for parquet files")
    parser.add_argument("--species", required=True, help="Species key")
    parser.add_argument("--min-tpm", type=float, default=0.0, help="Filter out records with TPM below this value")
    parser.add_argument("--validate-genes", action="store_true", help="Validate gene IDs against geneinfo.db")
    parser.add_argument("--dry-run", action="store_true", help="Run validation and parsing without writing output")
    args = parser.parse_args()

    converter = BulkRNARsemConverter(
        input_dir=Path(args.input_dir),
        metadata_file=Path(args.metadata),
        output_dir=Path(args.output_dir),
        species=args.species,
        min_tpm=args.min_tpm,
        validate_genes=args.validate_genes,
    )
    success = converter.run(dry_run=args.dry_run)
    return 0 if success else 1


if __name__ == "__main__":
    raise SystemExit(main())
