#!/usr/bin/env python3
"""
scATAC Gene Activity -> h5ad converter (aligned with the scRNA h5ad pipeline).

Highlights:
- Input: Signac/Seurat `.rds` / `.RData` containing GeneActivity results.
- Output: one h5ad per sample (X=log1p Gene Activity, CSC+gzip) + `*_metadata.json`.
- Keeps only GAS / obs / UMAP / var; optionally records fragment paths for future coverage features.
- No silent fallback: missing GeneActivity / clusters / UMAP will fail fast.
- Strict by default: only accepts the specified assay/slot log1p matrix; use `--allow-counts` to normalize+log1p from counts.
"""
import os
import sys
import argparse
import json
import time
import logging
from pathlib import Path
from typing import List, Optional, Tuple, Dict

import numpy as np
import pandas as pd
import scipy.sparse as sp
import anndata as ad
from pandas.api import types as pdt

# logging
SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))
try:
    from logging_utils import configure_script_logger, create_print_proxy
    _LOGGER = configure_script_logger(__name__)
    print = create_print_proxy(_LOGGER)
except Exception:
    logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
    _LOGGER = logging.getLogger(__name__)
    print = _LOGGER.info

# Prefer conda R
_conda_r_home = Path(sys.executable).resolve().parent.parent / "lib" / "R"
if "R_HOME" not in os.environ and _conda_r_home.exists():
    os.environ["R_HOME"] = str(_conda_r_home)

try:
    import rpy2  # noqa: F401
except ImportError:
    rpy2 = None

ro = None
pandas2ri = None
default_converter = None
localconverter = None
callbacks = None


def _require_rpy2() -> None:
    """Import rpy2 lazily so `--help` works without heavy dependencies."""

    global ro, pandas2ri, default_converter, localconverter, callbacks
    if ro is not None:
        return
    if rpy2 is None:
        raise RuntimeError("rpy2 is required. Install it first (e.g. `mamba install rpy2`).")

    import rpy2.robjects as _ro
    from rpy2.robjects import pandas2ri as _pandas2ri, default_converter as _default_converter
    from rpy2.robjects.conversion import localconverter as _localconverter
    from rpy2.rinterface_lib import callbacks as _callbacks

    ro = _ro
    pandas2ri = _pandas2ri
    default_converter = _default_converter
    localconverter = _localconverter
    callbacks = _callbacks


# ---- helpers -----------------------------------------------------

def _sanitize_obs(df: pd.DataFrame) -> pd.DataFrame:
    """Ensure obs columns are HDF5-writable."""
    from pandas import CategoricalDtype
    for col in df.columns:
        dt = df[col].dtype
        if pdt.is_extension_array_dtype(dt) or isinstance(dt, CategoricalDtype) or dt == object:
            if pdt.is_numeric_dtype(dt):
                df[col] = df[col].astype(float)
            else:
                df[col] = df[col].astype(str)
    return df


def _align_index(df: pd.DataFrame, names: List[str], label: str) -> pd.DataFrame:
    """Strictly align indexes; fail on content mismatch or duplicates."""
    tgt = list(map(str, names))
    src = df.index.astype(str).tolist()

    if set(tgt) != set(src):
        missing = set(tgt) - set(src)
        extra = set(src) - set(tgt)
        raise ValueError(
            f"{label} index mismatch: missing {len(missing)}, extra {len(extra)}"
        )
    if df.index.has_duplicates:
        dup = df.index[df.index.duplicated()].unique().tolist()
        raise ValueError(f"{label} index contains duplicates, e.g. {dup[:5]}")
    aligned = df.reindex(tgt)
    # If there are no columns (e.g. meta.features only has row names), allow it.
    if aligned.shape[1] > 0 and aligned.isnull().all(axis=1).any():
        bad = aligned.index[aligned.isnull().all(axis=1)][:5].tolist()
        raise ValueError(f"{label} alignment produced fully-missing rows, e.g. {bad}")
    return aligned


def _normalize_log1p(X: sp.csr_matrix, target_sum: float = 1e4) -> sp.csr_matrix:
    """Library-size normalize per cell and log1p (in-place)."""
    counts = np.asarray(X.sum(axis=1)).ravel()
    counts[counts == 0] = 1
    scale = target_sum / counts
    X.data *= np.repeat(scale, np.diff(X.indptr))
    np.log1p(X.data, out=X.data)
    return X


class ScATACH5ADConverter:
    def __init__(
        self,
        input_path: Path,
        output_root: Path,
        overwrite: bool = False,
        assay: str = "ACTIVITY",
        slot: str = "data",
        allow_counts: bool = False,
        cluster_col: Optional[str] = None,
        fragment_path_col: Optional[str] = None,
        fragment_path_default: Optional[str] = None,
        compression: str = "gzip",
        write_metadata_hash: bool = False,
        object_name: Optional[str] = None,
        umap_name: Optional[str] = None,
    ):
        self.input_path = input_path
        self.output_root = output_root
        self.overwrite = overwrite
        self.assay = assay
        self.slot = slot
        self.allow_counts = allow_counts
        self.cluster_col = cluster_col
        self.fragment_path_col = fragment_path_col
        self.fragment_path_default = fragment_path_default
        self.compression = compression
        self.write_metadata_hash = write_metadata_hash
        self.object_name = object_name
        self.umap_name = umap_name
        self.r_initialized = False

    # ---- R init & loaders ----
    def _init_r(self):
        if self.r_initialized:
            return
        _require_rpy2()
        callbacks.consolewrite_print = lambda x: None
        callbacks.consolewrite_warn = lambda x: None
        self.r = ro.r
        self.r('suppressMessages(library(Seurat))')
        self.r('suppressMessages(library(Signac))')
        self.r_initialized = True

    def _load_r_object(self, path: Path):
        self._init_r()
        r_path = str(path).replace('"', '\\"')
        cls = []
        if path.suffix.lower() == ".rds":
            self.r(f'obj <- readRDS("{r_path}")')
            cls = list(self.r('class(obj)'))
        else:
            self.r(f'loaded_objs <- load("{r_path}")')
            self.r('obj <- NULL')
            if self.object_name:
                self.r(f'if ("{self.object_name}" %in% loaded_objs) obj <- get("{self.object_name}")')
            else:
                self.r('if (length(loaded_objs) == 1) { obj <- get(loaded_objs[[1]]) }')
            cls = list(self.r('class(obj)')) if not self.r('is.null(obj)')[0] else []

        if "Seurat" not in cls:
            loaded = []
            try:
                loaded = list(self.r('if (exists("loaded_objs")) loaded_objs else character(0)'))
            except Exception:
                loaded = []
            msg = f"File {path} is not a Seurat/Signac object (class={cls}, loaded={loaded})"
            if path.suffix.lower() != ".rds":
                msg += "; use --object-name to select the object when a .RData contains multiple objects"
            raise ValueError(msg)
        # Cleanup temporary R objects to avoid cross-sample contamination.
        self.r('rm(list=setdiff(ls(), c("obj")))')
        return "obj"

    # ---- extraction helpers ----
    def _extract_matrix(self, r_expr: str, force_int: bool = False) -> Optional[Tuple[sp.csr_matrix, List[str], List[str]]]:
        """Robust extraction preferring dgCMatrix (R is genes×cells; convert to cells×genes CSR)."""
        self.r(f'suppressWarnings(tmp <- {r_expr})')
        if self.r('is.null(tmp)')[0]:
            return None
        try:
            self.r('suppressWarnings(if (!inherits(tmp, "dgCMatrix")) tmp <- as(tmp, "dgCMatrix"))')
            with localconverter(default_converter + pandas2ri.converter):
                data = np.array(self.r('tmp@x'), dtype=np.float32)
                indices = np.array(self.r('tmp@i'), dtype=np.int32)
                indptr = np.array(self.r('tmp@p'), dtype=np.int32)
                dims = list(self.r('dim(tmp)'))
                rows = list(self.r('rownames(tmp)'))
                cols = list(self.r('colnames(tmp)'))
            self.r('rm(tmp)')
            mat = sp.csc_matrix((data, indices, indptr), shape=(dims[0], dims[1])).T.tocsr()
            if force_int:
                mat.data = mat.data.astype(np.int32)
            return mat, rows, cols
        except Exception as e:
            print(f"   ⚠️ Matrix extraction failed: {e}")
            self.r('rm(tmp)')
            return None

    # ---- core logic ----
    def _ensure_assay_slot(self, obj_sym: str) -> str:
        assays = list(self.r(f'Assays({obj_sym})'))
        if not assays:
            raise ValueError("No assays found in the object")
        if self.assay not in assays:
            raise ValueError(f"Requested assay={self.assay} not found. Available: {assays}")
        slot_ok = bool(self.r(f'!is.null(tryCatch(GetAssayData({obj_sym}, assay=\"{self.assay}\", slot=\"{self.slot}\"), error=function(e) NULL))')[0])
        if not slot_ok:
            if self.allow_counts and self.slot != "counts":
                # Allow fallback to counts later, but first ensure counts exists.
                counts_ok = bool(self.r(f'!is.null(tryCatch(GetAssayData({obj_sym}, assay=\"{self.assay}\", slot=\"counts\"), error=function(e) NULL))')[0])
                if not counts_ok:
                    raise ValueError(f"assay={self.assay} missing slot={self.slot}, and counts is also missing; cannot continue")
            else:
                raise ValueError(
                    f"assay={self.assay} missing slot={self.slot}; run NormalizeData/log1p or adjust --slot/--allow-counts"
                )
        return self.assay

    def _get_cluster_col(self, obs: pd.DataFrame) -> Tuple[str, List[str]]:
        """
        Keep all groupable columns by default (string/categorical, or numeric with <=100 unique values).
        Returns (default_col, candidate_cols). If none found, inject a placeholder `cluster_all`.
        """
        if self.cluster_col:
            if self.cluster_col not in obs.columns:
                raise ValueError(f"Requested cluster column {self.cluster_col} not found. Candidates: {list(obs.columns)}")
            return self.cluster_col, [self.cluster_col]

        blocklist = {
            "n_genes", "n_counts", "total_counts", "n_counts_rna", "n_features_rna",
            "percent_mito", "scrublet_score", "doublet_score", "s_score", "g2m_score",
            "batch", "phase", "barcode", "cell_id",
            # Common ATAC QC columns
            "peak_region_fragments", "n_fragments", "ncount_atac", "tss_enrichment",
            "pct_reads_in_peaks", "blacklist_ratio", "passed_filters",
        }

        candidates: List[str] = []
        for col in obs.columns:
            series = obs[col]
            col_lower = col.lower()
            if col_lower in blocklist:
                continue
            if isinstance(series.dtype, pd.CategoricalDtype) or pd.api.types.is_object_dtype(series):
                candidates.append(col)
                continue
            if pd.api.types.is_integer_dtype(series) or pd.api.types.is_float_dtype(series):
                try:
                    n_unique = series.nunique(dropna=True)
                    if n_unique <= 100:
                        candidates.append(col)
                except Exception:
                    continue

        if not candidates:
            placeholder = "cluster_all"
            obs[placeholder] = "all"
            return placeholder, [placeholder]

        # Default to the first candidate (input column order).
        return candidates[0], candidates

    def _resolve_fragment_path(self, obs: pd.DataFrame) -> Optional[str]:
        if self.fragment_path_col and self.fragment_path_col in obs.columns:
            vals = obs[self.fragment_path_col].dropna().unique().tolist()
            if vals:
                return str(vals[0])
        if self.fragment_path_default:
            return self.fragment_path_default
        # Also try obj@misc as a fallback.
        try:
            frag = str(self.r('if (!is.null(obj@misc$fragments)) obj@misc$fragments else NA')[0])
            if frag != "NA":
                return frag
        except Exception:
            pass
        return None

    def _sha1_list(self, items: List[str]) -> str:
        import hashlib
        h = hashlib.sha1()
        for item in items:
            h.update(str(item).encode("utf-8"))
            h.update(b"\n")
        return h.hexdigest()

    def _validate_data_matrix(self, X: sp.csr_matrix):
        """Fail-fast guards for an already-normalized/log1p matrix."""
        if X.nnz == 0:
            raise ValueError("Data matrix is empty (nnz=0)")
        if X.data.size == 0:
            raise ValueError("Data matrix has no values")
        max_val = float(X.data.max())
        min_val = float(X.data.min())
        # Sample-based integer check to avoid scanning large matrices.
        sample_size = min(1000, X.data.size)
        sample_data = X.data[:sample_size]
        if sample_size > 0 and np.all(np.mod(sample_data, 1) == 0):
            if max_val > 50 or sample_size == X.data.size:
                raise ValueError(
                    "Detected integer-like values; likely raw counts. Run NormalizeData+log1p (or use --allow-counts)."
                )
        if max_val > 50:
            raise ValueError(f"Max value {max_val:.2f} is too large; likely not log-normalized (threshold: 50)")
        if min_val < 0:
            raise ValueError(
                f"Matrix contains negative values ({min_val:.2f}); likely scaled/residual. Provide log1p GeneActivity."
            )

    def _extract_obs(self, obj_sym: str) -> pd.DataFrame:
        with localconverter(default_converter + pandas2ri.converter):
            obs = pd.DataFrame(self.r(f'{obj_sym}@meta.data'))
        obs.index = obs.index.astype(str)
        return _sanitize_obs(obs)

    def _extract_var(self, assay: str) -> pd.DataFrame:
        # Note: meta.features may have 0 columns but valid row names; keep row names explicitly.
        # Seurat v5 compatibility: Assay5 doesn't have meta.features slot
        gene_names = []
        var = pd.DataFrame()

        # Try to get gene names first (works for both v3/v4 and v5)
        try:
            # First try Features() function (Seurat v5 compatible)
            gene_names = list(self.r(f'Features(obj[[\"{assay}\"]])'))
        except Exception:
            try:
                # Fallback: rownames of meta.features (Seurat v3/v4)
                gene_names = list(self.r(f'rownames(obj[[\"{assay}\"]]@meta.features)'))
            except Exception:
                gene_names = []

        # Try to get meta.features (only exists in Seurat v3/v4)
        try:
            with localconverter(default_converter + pandas2ri.converter):
                var = pd.DataFrame(self.r(f'obj[[\"{assay}\"]]@meta.features'))
        except Exception:
            # Seurat v5: meta.features slot doesn't exist, create empty DataFrame
            var = pd.DataFrame()

        # Force R-side row names to avoid index loss when there are 0 columns.
        if gene_names:
            var.index = pd.Index([str(g) for g in gene_names])
        else:
            var.index = var.index.astype(str)
        return var

    def _extract_all_umaps(self, obj_sym: str) -> Tuple[Dict[str, pd.DataFrame], str, List[str]]:
        """
        Extract all UMAP reductions and determine default.

        Returns:
            Tuple of (all_umaps_dict, default_umap_name, all_umap_names)
        """
        # Flexible UMAP detection: find any reduction containing "umap" (case-insensitive)
        reduction_list = list(self.r(f'names({obj_sym}@reductions)'))
        umap_reductions = [r for r in reduction_list if 'umap' in r.lower()]

        if not umap_reductions:
            raise ValueError(f"No UMAP reduction found. Available reductions: {reduction_list}")

        # Determine default UMAP
        if len(umap_reductions) == 1:
            default_umap = umap_reductions[0]
            print(f"   ℹ️ Using UMAP: {default_umap}")
        elif self.umap_name:
            if self.umap_name not in umap_reductions:
                raise ValueError(
                    f"--umap-name '{self.umap_name}' not found.\n"
                    f"   Available UMAP reductions: {umap_reductions}"
                )
            default_umap = self.umap_name
            print(f"   ℹ️ Using default UMAP: {default_umap} (from --umap-name)")
        else:
            raise ValueError(
                f"Multiple UMAP reductions found: {umap_reductions}\n"
                f"   Please specify one with --umap-name.\n"
                f"   Example: --umap-name {umap_reductions[0]}"
            )

        if len(umap_reductions) > 1:
            print(f"   ℹ️ Multiple UMAP reductions available: {umap_reductions}")

        # Extract all UMAPs
        all_umaps = {}
        for umap_name in umap_reductions:
            with localconverter(default_converter + pandas2ri.converter):
                umap = pd.DataFrame(self.r(f'Embeddings({obj_sym}, reduction=\"{umap_name}\")'))

            # Bind row names
            try:
                r_names = list(self.r(f'rownames(Embeddings({obj_sym}, reduction=\"{umap_name}\"))'))
                if len(r_names) == len(umap):
                    umap.index = pd.Index([str(x) for x in r_names], name="cell_id")
                else:
                    print(f"   ⚠️ {umap_name}: rowname length mismatch: {len(r_names)} vs {len(umap)}")
            except Exception as e:
                print(f"   ⚠️ {umap_name}: Failed to extract row names: {e}")

            if umap.shape[1] < 2:
                raise ValueError(f"{umap_name} has fewer than 2 dimensions")

            umap.index = umap.index.astype(str)
            umap = umap.iloc[:, :2]
            umap.columns = ["UMAP1", "UMAP2"]

            all_umaps[umap_name] = umap
            print(f"   ✅ Extracted {umap_name} ({umap.shape[0]} cells)")

        return all_umaps, default_umap, umap_reductions

    def convert_one(self, r_path: Path, index: int, total: int, study: Optional[str] = None) -> str:
        """Convert a single `.rds`/`.RData` file and return a status string."""
        print("=" * 60)
        obj_sym = self._load_r_object(r_path)
        sample = r_path.stem
        study_name = study or r_path.parent.name

        # Output path: if study is not provided, mirror the input directory structure under output_root.
        if study:
            out_dir = self.output_root / study / sample
        else:
            base_dir = Path(".")
            if self.input_path.is_dir():
                try:
                    base_dir = r_path.parent.relative_to(self.input_path)
                except Exception:
                    base_dir = Path(".")
            out_dir = self.output_root / base_dir / sample
        out_dir.mkdir(parents=True, exist_ok=True)
        h5ad_path = out_dir / f"{sample}.h5ad"
        meta_path = out_dir / f"{sample}_metadata.json"

        print(f"[{index}/{total}] 🧬 Converting {study_name}/{sample}")

        if h5ad_path.exists() and not self.overwrite:
            print(f"⏭️ Skip existing {h5ad_path} (use --overwrite to rebuild)")
            return "skipped"

        t0 = time.time()
        assay = self._ensure_assay_slot(obj_sym)
        print(f"   🧪 Using assay: {assay}, slot: {self.slot}")

        mat_data = self._extract_matrix(f'GetAssayData({obj_sym}, assay=\"{assay}\", slot=\"{self.slot}\")')
        use_counts = False
        if (mat_data is None or mat_data[0].nnz == 0) and self.slot != "counts" and self.allow_counts:
            print("   ⚠️ Slot data is empty; using counts and normalizing (because --allow-counts is set)")
            mat_data = self._extract_matrix(f'GetAssayData({obj_sym}, assay=\"{assay}\", slot=\"counts\")', force_int=True)
            use_counts = True
        if mat_data is None:
            raise ValueError("Failed to extract the GeneActivity matrix for the requested slot")

        X, genes, cells = mat_data
        print(f"   Matrix cells×genes: {X.shape[0]} x {X.shape[1]} (source={'counts' if use_counts else self.slot})")
        if use_counts:
            X = _normalize_log1p(X)
        else:
            if X.nnz == 0:
                raise ValueError(
                    "Extracted slot matrix has nnz=0; ensure GeneActivity exists or use --allow-counts"
                )
            self._validate_data_matrix(X)
            # Ensure float32
            X.data = X.data.astype(np.float32, copy=False)
        X = X.tocsc()

        obs = self._extract_obs(obj_sym)
        obs = _align_index(obs, cells, "obs")
        cluster_col, cluster_candidates = self._get_cluster_col(obs)
        # Convert candidate cluster columns to category to reduce size and speed up loading.
        for col in cluster_candidates:
            try:
                obs[col] = obs[col].astype("category")
            except Exception:
                pass

        var = self._extract_var(assay)
        # If meta.features is missing or row names are wrong, fall back to a minimal var (row names only).
        # KISS: keep var index aligned with matrix columns.
        # noqa: E265
        var_aligned = True
        try:
            var = _align_index(var, genes, "var")
        except Exception as e:
            print(f"   ⚠️ var alignment failed; using minimal var (row names only). Reason: {e}")
            var = pd.DataFrame(index=pd.Index([str(g) for g in genes], name="gene_id"))
            var_aligned = False

        # Extract all UMAP reductions
        all_umaps, default_umap, all_umap_names = self._extract_all_umaps(obj_sym)

        # Align all UMAPs to cells
        aligned_umaps = {}
        umap_fallbacks = {}
        for umap_name, umap in all_umaps.items():
            try:
                umap_aligned = _align_index(umap, cells, f"UMAP ({umap_name})")
                aligned_umaps[umap_name] = umap_aligned
                umap_fallbacks[f"umap_reindexed_{umap_name}"] = False
            except Exception as e:
                print(f"   ⚠️ {umap_name}: alignment failed; reindexing to cells order. Reason: {e}")
                # If lengths match, assume order matches and overwrite the index.
                if len(umap) == len(cells):
                    umap_aligned = umap.copy()
                    umap_aligned.index = pd.Index([str(c) for c in cells])
                    aligned_umaps[umap_name] = umap_aligned
                else:
                    # Reindex by cells; fail if still missing.
                    umap_aligned = umap.reindex(pd.Index([str(c) for c in cells]))
                    missing = umap_aligned.isnull().all(axis=1)
                    if missing.any():
                        bad = umap_aligned.index[missing].tolist()[:5]
                        raise ValueError(f"{umap_name}: still has missing rows after alignment, e.g. {bad}")
                    aligned_umaps[umap_name] = umap_aligned
                umap_fallbacks[f"umap_reindexed_{umap_name}"] = True

        # Process fragment files
        frag_path = self._resolve_fragment_path(obs)
        frag_effective = None
        has_frag = False
        if frag_path:
            src_frag = Path(frag_path).resolve()
            src_idx = src_frag.with_suffix(src_frag.suffix + ".tbi")
            if not src_idx.exists():
                # Common index file name: fragments.tsv.gz.tbi
                alt_idx = Path(str(src_frag) + ".tbi")
                if alt_idx.exists():
                    src_idx = alt_idx
            try:
                if src_frag.exists():
                    dst_frag = out_dir / src_frag.name
                    if dst_frag.exists() or dst_frag.is_symlink():
                        dst_frag.unlink()
                    dst_frag.symlink_to(src_frag)
                    if src_idx.exists():
                        dst_idx = out_dir / src_idx.name
                        if dst_idx.exists() or dst_idx.is_symlink():
                            dst_idx.unlink()
                        dst_idx.symlink_to(src_idx)
                    frag_effective = src_frag.name
                    has_frag = True
                    print(f"   🔗 Linked fragment: {dst_frag}")
                else:
                    print(f"   ⚠️ Fragment file does not exist, skipping: {src_frag}")
            except Exception as e:
                print(f"   ⚠️ Failed to link fragment; keeping original path: {e}")
                frag_effective = str(frag_path)
                has_frag = True

        # Build metadata
        meta = {
            "sample_name": sample,
            "study": study_name,
            "n_cells": int(X.shape[0]),
            "n_genes": int(X.shape[1]),
            "cluster_col": cluster_col,
            "cluster_candidates": cluster_candidates,
            "has_fragment": has_frag,
            "files": {"h5ad": h5ad_path.name},
            "source": r_path.name,
            "generated_at": time.strftime("%Y-%m-%d %H:%M:%S"),
            "embeddings": {
                "default": default_umap,
                "available": all_umap_names,
            },
            "assay_used": assay,
            "slot_used": "counts+log1p" if use_counts else self.slot,
            "allow_counts": use_counts,
            "fallbacks": {
                "var_aligned": var_aligned,
                **umap_fallbacks,
            },
        }

        # Build obsm: X_umap (default) + X_{original_name} for all
        obsm = {}
        obsm["X_umap"] = aligned_umaps[default_umap].values  # Standard key for frontend
        for umap_name, umap_data in aligned_umaps.items():
            obsm[f"X_{umap_name}"] = umap_data.values

        adata = ad.AnnData(X=X, obs=obs, var=var, obsm=obsm)
        adata.uns["fragment_path"] = frag_effective
        adata.uns["cluster_col"] = cluster_col
        adata.write_h5ad(h5ad_path, compression=self.compression)
        if has_frag:
            meta["fragment_path"] = frag_effective
        if self.write_metadata_hash:
            meta["hash"] = {
                "cells_sha1": self._sha1_list(cells),
                "genes_sha1": self._sha1_list(genes),
            }
        meta_path.write_text(json.dumps(meta, indent=2))

        print(f"   ✅ Finished {sample} in {time.time()-t0:.2f}s -> {h5ad_path}")
        return "success"

    def run(self, study: Optional[str] = None):
        targets: List[Path]
        if self.input_path.is_file():
            targets = [self.input_path]
        else:
            targets = sorted(
                list(self.input_path.rglob("*.rds")) + list(self.input_path.rglob("*.Rdata")),
                key=lambda p: str(p)
            )
        if not targets:
            raise FileNotFoundError(f"No .rds/.RData files found under {self.input_path}")

        print("=" * 60)
        print("🧬 scATAC GeneActivity -> h5ad")
        print(f"Input     : {self.input_path}")
        print(f"Output    : {self.output_root}")
        if study:
            print(f"Study     : {study}")
        print(f"Overwrite : {self.overwrite}")
        print(f"Assay     : {self.assay}")
        print(f"Slot      : {self.slot}")
        print(f"Allow cnt : {self.allow_counts}")
        print(f"Hash meta : {self.write_metadata_hash}")
        print("=" * 60)

        print("\n📋 Task list:")
        for i, f in enumerate(targets, 1):
            rel = f.relative_to(self.input_path) if f.is_relative_to(self.input_path) else f
            print(f"  {i}. {rel}")
        print(f"  Total: {len(targets)}\n")

        stats = {"success": 0, "skipped": 0, "failed": 0}
        total = len(targets)
        for idx, r_file in enumerate(targets, 1):
            try:
                status = self.convert_one(r_file, idx, total, study=study)
            except Exception as e:
                print(f"   ❌ Conversion failed: {e}")
                status = "failed"
            stats[status or "failed"] = stats.get(status or "failed", 0) + 1

        print("\n" + "=" * 60)
        print("🎉 Done")
        print(f"  Total   : {total}")
        print(f"  Success : {stats['success']}")
        print(f"  Skipped : {stats['skipped']}")
        print(f"  Failed  : {stats['failed']}")
        print("=" * 60)


# ---- CLI ----

def parse_args():
    ap = argparse.ArgumentParser(description="scATAC GeneActivity -> h5ad converter")
    ap.add_argument("--input", required=True, help="Input .rds/.RData file or a directory to scan recursively")
    ap.add_argument("--output-dir", required=True, help="Output root directory")
    ap.add_argument("--study", help="Optional: override study name for the output layout")
    ap.add_argument("--overwrite", action="store_true", help="Overwrite existing outputs")
    ap.add_argument("--assay", default="ACTIVITY", help="Assay name that stores GeneActivity (default: ACTIVITY)")
    ap.add_argument("--slot", default="data", help="Slot name for GeneActivity (default: data; counts fallback requires --allow-counts)")
    ap.add_argument("--allow-counts", action="store_true", help="If slot is empty, allow counts fallback and normalize+log1p (disabled by default)")
    ap.add_argument("--cluster-col", help="Cluster column name in obs (optional)")
    ap.add_argument("--fragment-path-col", help="obs column name for fragment path (optional)")
    ap.add_argument("--fragment-path-default", help="Default fragment path when missing (optional)")
    ap.add_argument("--compression", default="gzip", help="h5ad compression (default: gzip)")
    ap.add_argument("--metadata-hash", action="store_true", help="Write SHA1 hashes for cells/genes into metadata for verification")
    ap.add_argument("--object-name", help="When input is .RData with multiple objects, explicitly select the Seurat object name")
    ap.add_argument("--umap-name", help="When multiple UMAP reductions exist, specify which one to use as default (required if multiple found)")
    return ap.parse_args()


def main():
    args = parse_args()
    input_path = Path(args.input).expanduser()
    output_root = Path(args.output_dir).expanduser()
    conv = ScATACH5ADConverter(
        input_path=input_path,
        output_root=output_root,
        overwrite=args.overwrite,
        assay=args.assay,
        slot=args.slot,
        allow_counts=args.allow_counts,
        cluster_col=args.cluster_col,
        fragment_path_col=args.fragment_path_col,
        fragment_path_default=args.fragment_path_default,
        compression=args.compression,
        write_metadata_hash=args.metadata_hash,
        object_name=args.object_name,
        umap_name=args.umap_name,
    )
    conv.run(study=args.study)


if __name__ == "__main__":
    main()
