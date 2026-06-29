#!/usr/bin/env python3
"""
Direct Seurat (.Rdata/.rds) to .h5ad converter, MAGIC, and cache feather backfiller.

Features:
- Subcommands: `convert` (full conversion), `magic` (add MAGIC to existing h5ad), `cache` (export UMAP+obs feather).
- Reads Seurat v3/v4/v5 objects via rpy2.
- Produces a single .h5ad file per sample.
- Strict Mode: Only converts integer counts unless forced.
- Robust matrix extraction (dgCMatrix + summary fallback).
- True In-place normalization (no Scanpy dependency, memory efficient).
- Full metadata support (obs, var from meta.features).
- Enhanced UX: Task lists, progress tracking, and summary reports.
- Separation of concerns: `convert` writes the base h5ad; `magic` writes a companion `_magic.h5ad`.
"""

import os
import sys
import argparse
import gc
import logging
import warnings
import time
import json
from pathlib import Path
from typing import Optional, Dict, List, Tuple, Any, Union, Sequence

import numpy as np
import pandas as pd
import scipy.sparse as sp
import anndata as ad
from pandas.api import types as pdt
try:
    import pyarrow as pa
    import pyarrow.feather as feather
except ImportError:
    pa = None
    feather = None

# Configure logging
SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

try:
    from logging_utils import configure_script_logger, create_print_proxy
    _LOGGER = configure_script_logger(__name__)
    print = create_print_proxy(_LOGGER)
except ImportError:
    logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
    _LOGGER = logging.getLogger(__name__)
    print = _LOGGER.info

# Prefer conda-provided R
_conda_r_home = Path(sys.executable).resolve().parent.parent / "lib" / "R"
if "R_HOME" not in os.environ and _conda_r_home.exists():
    os.environ["R_HOME"] = str(_conda_r_home)

robjects = None
pandas2ri = None
default_converter = None
localconverter = None
callbacks = None


def _require_rpy2() -> None:
    """Import rpy2 lazily so `--help` works without heavy dependencies."""

    global robjects, pandas2ri, default_converter, localconverter, callbacks
    if robjects is not None:
        return
    try:
        import rpy2.robjects as _robjects
        from rpy2.robjects import pandas2ri as _pandas2ri, default_converter as _default_converter
        from rpy2.robjects.conversion import localconverter as _localconverter
        from rpy2.rinterface_lib import callbacks as _callbacks
    except ImportError as exc:
        raise RuntimeError("rpy2 is required. Install it first (e.g. `mamba install rpy2`).") from exc

    robjects = _robjects
    pandas2ri = _pandas2ri
    default_converter = _default_converter
    localconverter = _localconverter
    callbacks = _callbacks


# --- Helpers ---

def _normalize_log1p(X: sp.csr_matrix, target_sum: float = 1e4) -> sp.csr_matrix:
    """Truly in-place library size normalization and log1p for CSR matrix."""
    # Calculate library size (sum per row)
    counts = np.array(X.sum(axis=1)).ravel()
    counts[counts == 0] = 1 # Avoid division by zero
    scaling = target_sum / counts
    
    # Apply scaling to data array
    # X.indptr defines the boundaries of rows in X.data
    # We repeat the scaling factor for each row to match the number of elements in that row
    X.data *= np.repeat(scaling, np.diff(X.indptr))
    
    # In-place log1p
    np.log1p(X.data, out=X.data)
    
    return X

def _sanitize_obs(df: pd.DataFrame) -> pd.DataFrame:
    """Ensure obs columns are writable to HDF5."""
    from pandas import CategoricalDtype
    for col in df.columns:
        dt = df[col].dtype
        if pdt.is_extension_array_dtype(dt) or isinstance(dt, CategoricalDtype) or dt == object:
            if pdt.is_numeric_dtype(dt):
                df[col] = df[col].astype(float)
            else:
                df[col] = df[col].astype(str)
    return df

def _align_index(df: pd.DataFrame, names: Sequence[str], label: str) -> pd.DataFrame:
    """Align dataframe index to provided names list. Raises ValueError on misalignment."""
    original_index_str = df.index.astype(str).tolist()
    target_index_str = list(map(str, names))

    # Strict check for exact match of content
    if set(original_index_str) != set(target_index_str):
        missing_in_original = set(target_index_str) - set(original_index_str)
        extra_in_original = set(original_index_str) - set(target_index_str)
        error_msg = f"Index content mismatch for {label}. "
        if missing_in_original:
            error_msg += f"Target has elements not in original ({len(missing_in_original)} missing; e.g., {list(missing_in_original)[:5]}). "
        if extra_in_original:
            error_msg += f"Original has elements not in target ({len(extra_in_original)} extra; e.g., {list(extra_in_original)[:5]})."
        raise ValueError(error_msg)
    
    if df.index.has_duplicates:
        duplicate_labels = df.index[df.index.duplicated()].unique().tolist()
        raise ValueError(f"Original {label} index contains duplicates ({len(duplicate_labels)} duplicates; e.g., {duplicate_labels[:5]}).")

    # Reindex
    aligned_df = df.reindex(target_index_str)

    # Allow sporadic missing values but forbid fully-missing rows (alignment bug).
    missing_rows = aligned_df.isnull().all(axis=1)
    if missing_rows.any():
        missing_labels = aligned_df.index[missing_rows].tolist()
        raise ValueError(
            f"Reindexing {label} introduced fully-missing rows ({len(missing_labels)}); "
            f"e.g., {missing_labels[:5]}"
        )
    
    if aligned_df.index.has_duplicates:
        duplicate_labels = aligned_df.index[aligned_df.index.duplicated()].unique().tolist()
        raise ValueError(f"Aligned {label} index contains duplicates after reindexing ({len(duplicate_labels)} duplicates; e.g., {duplicate_labels[:5]}).")

    return aligned_df


class ScRNAH5ADConverter:
    def __init__(
        self,
        input_dir: Path,
        output_dir: Path = None, # Optional for magic-only mode
        magic_params: Optional[Dict] = None,
        dry_run: bool = False,
        overwrite: bool = False,
        study_filter: Optional[str] = None,
        force_norm_from: Optional[str] = None,
        magic_memory_budget: float = 50.0,
        magic_trim_threshold: float = 1e-2,
        compression: str = "gzip",
        magic_compression: str = "gzip",
    ):
        self.input_dir = input_dir.resolve()
        self.output_dir = output_dir.resolve() if output_dir else None
        self.magic_params = magic_params or {}
        self.dry_run = dry_run
        self.overwrite = overwrite
        self.study_filter = study_filter
        self.force_norm_from = force_norm_from
        
        # MAGIC defaults
        self.magic_default_chunk_size = 5000
        self.magic_trim_threshold = magic_trim_threshold
        self.magic_memory_budget_gb = magic_memory_budget

        self.compression = compression
        self.magic_compression = magic_compression
        
        # Initialize R only if needed (convert mode)
        self.r_initialized = False

    def _init_r(self):
        """Initialize R instance and load Seurat."""
        if self.r_initialized:
            return
        _require_rpy2()
        try:
            # Strong suppression of R output
            callbacks.consolewrite_print = lambda x: None
            callbacks.consolewrite_warn = lambda x: None
            self.r = robjects.r
            # No initial print, keep clean
            self.r('suppressMessages(library(Seurat))')
            self.r('suppressMessages(library(methods))')
            self.r_initialized = True
        except Exception as e:
            raise RuntimeError(f"Failed to initialize R/Seurat: {e}") from e

    # ---
    # Scanning ---

    def scan_r_inputs(self) -> List[Tuple[Path, str, Path]]:
        """Recursively find .Rdata/.rds files."""
        files = []
        if not self.input_dir.exists():
            return []
        for p in self.input_dir.rglob("*"):
            if p.is_file() and p.suffix.lower() in (".rdata", ".rds"):
                if self.study_filter and self.study_filter not in str(p): continue
                try:
                    rel = p.relative_to(self.input_dir)
                    # Use parent folder as study name
                    study_name = rel.parent
                    if str(study_name) == ".": study_name = Path("root")
                    files.append((study_name, p.stem, p.resolve()))
                except: continue
        files.sort()
        return files

    def scan_h5ad_inputs(self) -> List[Path]:
        """Recursively find .h5ad files for MAGIC backfill."""
        files = []
        if not self.input_dir.exists():
            return []
        for p in self.input_dir.rglob("*.h5ad"):
            if self.study_filter and self.study_filter not in str(p): continue
            name_lower = p.name.lower()
            if name_lower.endswith("_magic.h5ad") or name_lower.endswith(".tmp"):  # skip existing MAGIC outputs / temp files
                continue
            stem = p.stem
            meta_path = p.with_name(f"{stem}_metadata.json")
            magic_path = p.with_name(f"{stem}_magic.h5ad")
            # Prefer cheap skips before scheduling
            if magic_path.exists():
                continue
            if meta_path.exists():
                try:
                    meta = json.loads(meta_path.read_text())
                    if "magic" in (meta.get("files") or {}):
                        continue
                    versions = meta.get("available_data_versions") or meta.get("versions") or {}
                    if isinstance(versions, dict) and "magic" in versions:
                        continue
                    if isinstance(versions, list) and "magic" in versions:
                        continue
                except Exception:
                    pass
            files.append(p.resolve())
        files.sort()
        return files

    def scan_cache_inputs(self) -> List[Path]:
        """Recursively find .h5ad files for cache feather export."""
        files = []
        if not self.input_dir.exists():
            return []
        for p in self.input_dir.rglob("*.h5ad"):
            if self.study_filter and self.study_filter not in str(p):
                continue
            name_lower = p.name.lower()
            if name_lower.endswith(".tmp") or name_lower.endswith("_magic.h5ad"):
                continue
            stem = p.stem
            meta_path = p.with_name(f"{stem}_metadata.json")
            cache_path = p.with_name(f"{stem}_cache.feather")

            if not self.overwrite:
                # Skip if cache already present on disk
                if cache_path.exists():
                    continue
                if meta_path.exists():
                    try:
                        meta = json.loads(meta_path.read_text())
                        files_meta = meta.get("files") or {}
                        cache_file = files_meta.get("cache_umap_obs")
                        if cache_file and (p.parent / cache_file).exists():
                            continue
                        cache_meta = meta.get("cache_meta") or {}
                        cache_file_meta = cache_meta.get("file")
                        if cache_file_meta and (p.parent / cache_file_meta).exists():
                            continue
                    except Exception:
                        pass
            files.append(p.resolve())
        files.sort()
        return files

    # ---
    # R Operations ---

    def _load_r_object(self, file_path: Path) -> str:
        """Load R file, return object name."""
        path_str = str(file_path).replace("\\", "/")
        try:
            if file_path.suffix.lower() == ".rds":
                print(f"   📖 Loading RDS: {file_path.name}")
                self.r.assign("seurat_obj", self.r(f'readRDS("{path_str}")'))
                return "seurat_obj"
            else:
                print(f"   📖 Loading RData: {file_path.name}")
                self.r(f'loaded_objs <- load("{path_str}")')
                objs = list(self.r('loaded_objs'))
                # Heuristic: find Seurat object
                for name in objs:
                     if self.r(f'inherits({name}, "Seurat")')[0]: return name
                     # Check inside list
                     if self.r(f'is.list({name})')[0]:
                         for sub in list(self.r(f'names({name})')):
                             if self.r(f'inherits({name}[["{sub}"]], "Seurat")')[0]:
                                 self.r(f'{name}_sub <- {name}[["{sub}"]]')
                                 return f"{name}_sub"
                raise ValueError(f"No Seurat object in {file_path.name}")
        except Exception as e:
            raise RuntimeError(f"Load failed: {e}")

    def _extract_matrix_robust(self, r_expr: str, force_int32: bool = False) -> Optional[Tuple[sp.csr_matrix, List[str], List[str]]]:
        """Extract matrix with dgCMatrix preference and summary() fallback."""
        # Suppress warnings during extraction
        self.r(f'suppressWarnings(tmp <- {r_expr})')
        if self.r('is.null(tmp)')[0]: return None
        
        # Try direct dgCMatrix extraction
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
            
            if len(data) == 0:
                mat = sp.csr_matrix((dims[1], dims[0]), dtype=np.float32)
            else:
                mat = sp.csc_matrix((data, indices, indptr), shape=(dims[0], dims[1])).T.tocsr()
            
            if force_int32: mat.data = mat.data.astype(np.int32)
            return mat, rows, cols
        except Exception:
            print(f"   ⚠️  Direct extraction failed. Trying summary()...")
            # Fallback to summary() for huge/weird matrices
            try:
                self.r('s <- summary(tmp)')
                with localconverter(default_converter + pandas2ri.converter):
                    i = np.array(self.r('s$i'), dtype=np.int32) - 1 # R 1-based
                    j = np.array(self.r('s$j'), dtype=np.int32) - 1
                    x = np.array(self.r('s$x'), dtype=np.float32)
                    dims = list(self.r('dim(tmp)'))
                    rows = list(self.r('rownames(tmp)'))
                    cols = list(self.r('colnames(tmp)'))
                
                # summary returns COO like (i=row, j=col) for R matrix (Genes x Cells)
                # We want Cells x Genes (transpose)
                # So R_i (gene) becomes Col, R_j (cell) becomes Row
                mat = sp.coo_matrix((x, (j, i)), shape=(dims[1], dims[0])).tocsr()
                if force_int32: mat.data = mat.data.astype(np.int32)
                self.r('rm(tmp); rm(s)')
                return mat, rows, cols
            except Exception as e2:
                print(f"   ❌ Summary extraction also failed: {e2}")
                self.r('rm(tmp)')
                return None

    def _get_metadata(self, obj_name: str, cell_names: List[str]) -> pd.DataFrame:
        """Extract metadata and align."""
        with localconverter(default_converter + pandas2ri.converter):
            try:
                df = robjects.conversion.rpy2py(self.r(f'{obj_name}@meta.data'))
            except:
                df = pd.DataFrame(index=cell_names)
        
        if not isinstance(df, pd.DataFrame):
             df = pd.DataFrame(index=cell_names)

        df = _align_index(df, cell_names, "obs")
        df = _sanitize_obs(df)
        df.index.name = "cell_id"
        return df

    def _get_gene_meta(self, obj_name: str, assay: str, gene_names: List[str]) -> pd.DataFrame:
        """Extract gene metadata from meta.features."""
        try:
            has_meta = bool(self.r(f'"meta.features" %in% slotNames({obj_name}[["{assay}"]])')[0])
        except Exception:
            has_meta = False
        if not has_meta:
            return pd.DataFrame(index=gene_names)

        with localconverter(default_converter + pandas2ri.converter):
            try:
                df = robjects.conversion.rpy2py(self.r(f'{obj_name}[["{assay}"]]@meta.features'))
            except:
                df = None
        
        if not isinstance(df, pd.DataFrame) or df.empty:
             return pd.DataFrame(index=gene_names)
        
        return _sanitize_obs(_align_index(df, gene_names, "var"))

    def _probe_assay(self, obj_name: str, assay: str, has_layerdata: bool) -> Dict:
        """Check assay slots for counts/data presence and quality."""
        res = {'counts': None, 'data': None}
        
        # Helper to probe a specific slot/layer
        def _check_slot(is_counts: bool):
            slot_type = "counts" if is_counts else "data"
            try:
                # Construct command
                if has_layerdata and self.r(f'inherits({obj_name}[["{assay}"]], "Assay5")')[0]:
                    cmd = f'LayerData({obj_name}, assay="{assay}", layer="{slot_type}")'
                else:
                    cmd = f'GetAssayData({obj_name}, assay="{assay}", slot="{slot_type}")'
                
                # Assign to temp var
                self.r(f'p_tmp <- {cmd}')
                
                # Check for NULL
                if self.r('is.null(p_tmp)')[0]: return None
                
                # Force dgCMatrix
                self.r('if (!inherits(p_tmp, "dgCMatrix")) p_tmp <- as(p_tmp, "dgCMatrix")')
                
                # Get samples (take head instead of random sample to be safe and fast)
                # If empty matrix, x is empty
                self.r('x_vals <- p_tmp@x')
                self.r('if (length(x_vals) > 1000) x_vals <- head(x_vals, 1000)')
                
                arr = np.array(self.r('x_vals'), dtype=np.float32)
                self.r('rm(p_tmp); rm(x_vals)')
                
                if arr.size > 0:
                    return {
                        'non_neg': bool((arr >= 0).all()),
                        'int_ratio': float((np.isclose(arr, np.round(arr)).sum() / arr.size))
                    }
                else:
                    # Empty matrix is valid (non-negative, integer)
                    return {'non_neg': True, 'int_ratio': 1.0}
            except Exception:
                return None

        res['counts'] = _check_slot(True)
        res['data'] = _check_slot(False)
        
        return res

    def _select_best_assay(self, obj_name: str) -> Tuple[str, bool]:
        """
        Select best assay (Strict Mode).
        Returns (assay, is_counts).
        """
        assays = list(self.r(f'names({obj_name}@assays)'))
        has_layerdata = bool(self.r('exists("LayerData")')[0])
        
        # Check force
        if self.force_norm_from:
            if self.force_norm_from in assays:
                p = self._probe_assay(obj_name, self.force_norm_from, has_layerdata)
                # Try to force as Counts first (if non-negative)
                if p['counts'] and p['counts']['non_neg']:
                    print(f"   ⚠️  Forcing assay '{self.force_norm_from}' as Counts (int_ratio={p['counts']['int_ratio']:.2f}).")
                    return self.force_norm_from, True
                print(f"   ❌ Forced assay '{self.force_norm_from}' has no valid non-negative counts slot.")
            else:
                print(f"   ❌ Forced assay '{self.force_norm_from}' not found.")

        # Strict: only consider RNA counts
        if "RNA" not in assays:
            raise ValueError("RNA assay not found; only RNA counts are supported.")

        p = self._probe_assay(obj_name, "RNA", has_layerdata)
        c = p['counts']
        if c and c['non_neg'] and c['int_ratio'] >= 0.98:
            print(f"   ✅ Selected 'RNA' (Counts, int_ratio={c['int_ratio']:.2f})")
            return "RNA", True

        if c:
            raise ValueError(f"RNA counts not integer-like (int_ratio={c['int_ratio']:.2f}).")

        raise ValueError("No usable RNA counts found.")

    # ---
    # MAGIC & Processing ---

    def _run_magic_memsafe(self, adata: ad.AnnData) -> Optional[sp.csr_matrix]:
        """Run MAGIC with strict memory checking."""
        n_cells, n_genes = adata.shape
        max_cells = self.magic_params.get('max_cells', 50000)
        
        if n_cells > max_cells:
            print(f"   ⏭️  Skip MAGIC: {n_cells} cells > {max_cells} limit")
            return None

        # Est. Peak: Dense Input (float32) + MAGIC output (float64)
        peak_gb = (n_cells * n_genes * (4 + 8)) / (1024**3)
        
        if peak_gb > self.magic_memory_budget_gb:
            print(f"   ❌ Skip MAGIC: Est. memory {peak_gb:.1f} GB exceeds budget {self.magic_memory_budget_gb} GB")
            return None
        
        print(f"   🪄 Running MAGIC (Est. {peak_gb:.1f} GB)...")

        try:
            import magic
            # Use chunks for densification if large
            if sp.issparse(adata.X):
                if (n_cells*n_genes*4/1e9) > 2.0:
                    dense = np.zeros((n_cells, n_genes), dtype=np.float32)
                    chunk = 5000
                    for i in range(0, n_cells, chunk):
                        e = min(i+chunk, n_cells)
                        dense[i:e, :] = adata.X[i:e].toarray()
                else:
                    dense = adata.X.toarray()
            else:
                dense = adata.X

            op = magic.MAGIC(solver='approximate', n_jobs=1, verbose=False)
            res = op.fit_transform(dense)
            del dense; gc.collect()
            
            res_csr = sp.csr_matrix(res, dtype=np.float32)
            
            if self.magic_trim_threshold > 0:
                mask = res_csr.data < self.magic_trim_threshold
                if mask.any():
                    res_csr.data[mask] = 0
                    res_csr.eliminate_zeros()
            
            return res_csr
        except Exception as e:
            print(f"   ❌ MAGIC failed: {e}")
            return None

    def _check_h5ad_integrity(self, path: Path) -> bool:
        """Quickly verify if H5AD can be opened."""
        try:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                ad.read_h5ad(path, backed='r').file.close()
            return True
        except Exception:
            return False

    def _write_h5ad_atomic(self, adata: ad.AnnData, final_path: Path, compression: str = "gzip"):
        """Write to temp file then rename to ensure atomicity."""
        temp_path = final_path.with_suffix(final_path.suffix + ".tmp")
        try:
            print(f"   💾 Writing {final_path.name}...")
            adata.write_h5ad(temp_path, compression=compression)
            temp_path.rename(final_path)
        except Exception as e:
            print(f"   ❌ Write failed: {e}")
            if temp_path.exists():
                temp_path.unlink()
            raise e

    @staticmethod
    def _hash_obs_names(names: Sequence[str]) -> str:
        import hashlib
        return hashlib.sha1("\n".join([str(n) for n in names]).encode("utf-8")).hexdigest()

    def _select_embedding_for_cache(self, adata: ad.AnnData) -> Optional[str]:
        """Pick an embedding key for cache export; prefer UMAP, fallback to TSNE/others."""
        priorities = ["X_umap", "X_tsne"]
        for key in priorities:
            if key in adata.obsm:
                return key
        # fallback: first 2D embedding
        for key, val in adata.obsm.items():
            if hasattr(val, "shape") and val.shape[1] >= 2:
                return key
        return None

    def _export_cache_feather(
        self,
        adata: ad.AnnData,
        out_dir: Path,
        sample_name: str,
        obs_cols: Optional[List[str]] = None,
        embed_key: Optional[str] = None,
    ) -> Tuple[Optional[str], Dict]:
        """
        Export a single Feather containing cell_id + UMAP coords + selected obs columns.
        Returns (filename or None, meta dict).
        """
        if pa is None or feather is None:
            return None, {}
        obs_names = [str(x) for x in adata.obs_names]
        obs_hash = self._hash_obs_names(obs_names)

        cols = [pa.array(obs_names, type=pa.string())]
        col_names = ["cell_id"]

        umap_cols = []
        embed_used = embed_key or ("X_umap" if "X_umap" in adata.obsm else None)
        if embed_used and embed_used in adata.obsm:
            emb = adata.obsm[embed_used]
            if isinstance(emb, sp.spmatrix):
                emb = emb.toarray()
            base = embed_used.replace("X_", "") if embed_used.startswith("X_") else embed_used
            base = base.lower() or "embed"
            for i in range(emb.shape[1]):
                col_name = f"{base}_{i+1}"
                col_names.append(col_name)
                cols.append(pa.array(np.asarray(emb[:, i]).astype(np.float32)))
                umap_cols.append(col_name)

        if obs_cols is None:
            obs_cols = list(adata.obs.columns)
        for col in obs_cols:
            series = adata.obs[col]
            # categories -> string while keeping the original values
            if str(series.dtype) == "category":
                series = series.astype(str)
            col_names.append(col)
            cols.append(pa.array(series.tolist()))

        table = pa.table(cols, names=col_names)
        fpath = out_dir / f"{sample_name}_cache.feather"
        feather.write_feather(table, fpath, compression="zstd")

        meta = {
            "file": fpath.name,
            "columns": col_names,
            "umap_columns": umap_cols,
            "obs_columns": obs_cols,
            "obs_hash": obs_hash,
            "embedding_key": embed_used,
            "format": "feather",
            "generated_at": time.strftime("%Y-%m-%d %H:%M:%S"),
            "version": 1,
        }
        return fpath.name, meta

    def process_conversion(self, study_dir: Path, sample_name: str, r_file: Path, index: int, total: int) -> str:
        """
        Full conversion pipeline.
        Returns status: "success", "skipped", "failed"
        """
        self._init_r()
        print("=" * 60)
        print(f"[{index}/{total}] 🚀 Processing: {study_dir}/{sample_name}")
        
        out_dir = self.output_dir / study_dir / sample_name
        h5ad_path = out_dir / f"{sample_name}.h5ad"
        json_path = out_dir / f"{sample_name}_metadata.json"
        
        if h5ad_path.exists() and not self.overwrite:
            if self._check_h5ad_integrity(h5ad_path):
                print(f"⏭️  Skip: Output exists and is valid.")
                return "skipped"
            else:
                print(f"⚠️  Output exists but corrupted. Reprocessing.")

        try:
            start_t = time.time()
            obj = self._load_r_object(r_file)
            assay, is_counts = self._select_best_assay(obj)
            
            # Extract
            raw_csr, raw_genes, raw_cells = None, [], []
            norm_csr, norm_genes, norm_cells = None, [], []
            
            has_ld = bool(self.r('exists("LayerData")')[0])
            is_v5 = bool(self.r(f'inherits({obj}[["{assay}"]], "Assay5")')[0])

            # Try get Counts (if selected assay is counts)
            if is_counts:
                cmd = f'LayerData({obj}, assay="{assay}", layer="counts")' if (is_v5 and has_ld) else f'GetAssayData({obj}, assay="{assay}", slot="counts")'
                res = self._extract_matrix_robust(cmd, force_int32=True)
                if res: raw_csr, raw_genes, raw_cells = res
            else:
                # Forced data path (rare; e.g., --force-norm-from)
                cmd = f'LayerData({obj}, assay="{assay}", layer="data")' if (is_v5 and has_ld) else f'GetAssayData({obj}, assay="{assay}", slot="data")'
                res = self._extract_matrix_robust(cmd, force_int32=False)
                if res: norm_csr, norm_genes, norm_cells = res
            
            # Align
            X = None
            layers = {}
            
            if not is_counts: # Data mode (Force)
                if norm_csr is None: raise ValueError("No data matrix found in forced assay.")
                X = norm_csr
                final_cells, final_genes = norm_cells, norm_genes
            else: # Counts mode
                if raw_csr is None: raise ValueError("No counts matrix found.")
                X = raw_csr.astype(np.float32) # Placeholder, will norm in-place
                final_cells, final_genes = raw_cells, raw_genes
                layers['raw'] = raw_csr
            
            # Metadata
            obs = self._get_metadata(obj, final_cells)
            var = self._get_gene_meta(obj, assay, final_genes)
            
            adata = ad.AnnData(X=X, obs=obs, var=var)
            for k,v in layers.items(): adata.layers[k] = v
            
            # Normalize if X was from raw
            if is_counts:
                print("   Computing Log1p...")
                _normalize_log1p(adata.X)

            # Store matrices as CSC for fast gene slicing
            if sp.issparse(adata.X):
                adata.X = adata.X.tocsc()
            for k in list(adata.layers.keys()):
                if sp.issparse(adata.layers[k]):
                    adata.layers[k] = adata.layers[k].tocsc()

            # Embeddings
            reductions = list(self.r(f'names({obj}@reductions)'))
            mapping = {'pca': 'X_pca', 'umap': 'X_umap', 'tsne': 'X_tsne'}
            for red in reductions:
                key = mapping.get(red.lower(), f"X_{red.lower()}")
                try:
                    # Capture output to suppress potential R spam
                    self.r(f'suppressWarnings(coords <- Embeddings({obj}, reduction="{red}"))')
                    coords = np.array(self.r('coords'), dtype=np.float32)
                    r_rows = list(self.r(f'rownames({obj}@reductions${red}@cell.embeddings)'))
                    if len(r_rows) != coords.shape[0]: continue
                    
                    emb_df = pd.DataFrame(coords, index=r_rows)
                    emb_df = _align_index(emb_df, final_cells, key)
                    adata.obsm[key] = emb_df.values
                except: continue

            # Metadata
            vers = {'norm': {'type': 'normalized', 'description': 'Log1p normalized'}}
            if 'raw' in adata.layers: vers['raw'] = {'type': 'raw_counts', 'description': 'Raw counts'}
            adata.uns['available_versions'] = vers
            adata.uns['default_version'] = 'norm'
            
            # Atomic Write
            if not self.dry_run:
                out_dir.mkdir(parents=True, exist_ok=True)
                self._write_h5ad_atomic(adata, h5ad_path, compression=self.compression)

                # Metadata.json
                meta = {
                    "sample_name": sample_name,
                    "study": str(study_dir),
                    "n_cells": adata.n_obs,
                    "n_genes": adata.n_vars,
                    "versions": list(vers.keys()),
                    "available_data_versions": vers,
                    "default_data_version": "norm",
                    "cell_metadata_columns": list(adata.obs.columns),
                    "embeddings": list(adata.obsm.keys()),
                    "embedding_files": {},
                    "embedding_meta": {},
                    "files": {"h5ad": h5ad_path.name},
                }
                with open(out_dir / f"{sample_name}_metadata.json", 'w') as f:
                    json.dump(meta, f, indent=2)
            
            print(f"   ✅ Done in {time.time()-start_t:.1f}s")
            del adata; gc.collect(); self.r('rm(list=ls()); gc()')
            return "success"
            
        except Exception as e:
            print(f"   ❌ Failed: {e}")
            # import traceback; traceback.print_exc()
            try: self.r('rm(list=ls()); gc()')
            except: pass
            return "failed"

    def process_magic_backfill(self, h5ad_path: Path, index: int, total: int) -> str:
        """Magic-only update."""
        print("=" * 60)
        print(f"[{index}/{total}] ✨ Processing: {h5ad_path.stem}")
        start_t = time.time()
        try:
            if not self._check_h5ad_integrity(h5ad_path):
                print(f"   ❌ Corrupted H5AD. Skipping.")
                return "failed"

            magic_path = h5ad_path.with_name(f"{h5ad_path.stem}_magic.h5ad")
            meta_path = h5ad_path.with_name(f"{h5ad_path.stem}_metadata.json")

            # Cheap skip: companion already exists
            if magic_path.exists():
                print(f"   ⏭️  Skip: {magic_path.name} already exists.")
                return "skipped"

            # Skip if metadata already declares magic (avoid double bookkeeping)
            if meta_path.exists():
                try:
                    meta = json.load(open(meta_path))
                    versions = meta.get("available_data_versions") or meta.get("versions") or {}
                    files = meta.get("files") or {}
                    if (isinstance(versions, dict) and "magic" in versions) or \
                       (isinstance(versions, list) and "magic" in versions) or \
                       ("magic_h5ad" in files):
                        print("   ⏭️  Skip: metadata already has MAGIC.")
                        return "skipped"
                except Exception:
                    pass

            adata = ad.read_h5ad(h5ad_path)

            # Always recompute MAGIC to ensure clean separation
            magic_matrix = self._run_magic_memsafe(adata)
            if magic_matrix is None:
                return "skipped"  # Skipped due to constraints or user limits

            if sp.issparse(magic_matrix):
                magic_csc = magic_matrix.tocsc().astype(np.float32)
            else:
                magic_csc = sp.csc_matrix(magic_matrix, dtype=np.float32)

            magic_adata = ad.AnnData(
                X=magic_csc,
                obs=adata.obs.copy(),
                var=adata.var.copy(),
                obsm=adata.obsm.copy(),
                uns={"available_versions": {"magic": {"type": "magic_imputed", "description": "MAGIC smoothed"}},
                     "default_version": "magic"},
            )

            if not self.dry_run:
                self._write_h5ad_atomic(magic_adata, magic_path, compression=self.magic_compression)

                # Update metadata (add magic companion)
                meta = {}
                if meta_path.exists():
                    try:
                        meta = json.load(open(meta_path))
                    except Exception:
                        meta = {}

                # Normalize structures
                meta.setdefault("files", {})
                meta["files"]["magic_h5ad"] = magic_path.name

                avail_versions = meta.get("available_data_versions")
                if not isinstance(avail_versions, dict):
                    avail_versions = {}
                avail_versions["magic"] = {"type": "magic_imputed", "description": "MAGIC smoothed"}
                meta["available_data_versions"] = avail_versions

                if "versions" in meta:
                    if isinstance(meta["versions"], list):
                        if "magic" not in meta["versions"]:
                            meta["versions"].append("magic")
                    elif isinstance(meta["versions"], dict):
                        meta["versions"]["magic"] = {"type": "magic_imputed", "description": "MAGIC smoothed"}

                meta.setdefault("default_data_version", "norm")

                with open(meta_path, 'w') as f:
                    json.dump(meta, f, indent=2)

            print(f"   ✅ MAGIC saved separately in {time.time()-start_t:.1f}s.")
            del adata; gc.collect()
            return "success"

        except Exception as e:
            print(f"   ❌ Backfill failed: {e}")
            return "failed"

    def process_cache_export(self, h5ad_path: Path, index: int, total: int) -> str:
        """Export UMAP+obs feather for an existing h5ad."""
        print("=" * 60)
        print(f"[{index}/{total}] 📦 Cache: {h5ad_path.stem}")
        try:
            meta_path = h5ad_path.with_name(f"{h5ad_path.stem}_metadata.json")
            cache_path = h5ad_path.with_name(f"{h5ad_path.stem}_cache.feather")

            if not self._check_h5ad_integrity(h5ad_path):
                print("   ❌ Corrupted H5AD. Skipping.")
                return "failed"

            if cache_path.exists() and not self.overwrite:
                print("   ⏭️  Skip: cache already exists.")
                return "skipped"

            meta = {}
            if meta_path.exists():
                try:
                    meta = json.loads(meta_path.read_text())
                    files_meta = meta.get("files") or {}
                    cache_meta = meta.get("cache_meta") or {}
                    if not self.overwrite:
                        cache_file = files_meta.get("cache_umap_obs")
                        if cache_file and (h5ad_path.parent / cache_file).exists():
                            print("   ⏭️  Skip: metadata already has cache entry.")
                            return "skipped"
                        cache_file_meta = cache_meta.get("file")
                        if cache_file_meta and (h5ad_path.parent / cache_file_meta).exists():
                            print("   ⏭️  Skip: cache file recorded in metadata exists.")
                            return "skipped"
                except Exception:
                    meta = {}

            adata = ad.read_h5ad(h5ad_path)
            embed_key = self._select_embedding_for_cache(adata)
            if embed_key is None:
                print("   ⚠️  No embedding found (UMAP/TSNE). Exporting obs only.")

            cache_file, cache_meta = self._export_cache_feather(
                adata, h5ad_path.parent, h5ad_path.stem, embed_key=embed_key
            )
            if cache_file is None:
                print("   ⚠️ Feather dependencies missing (pyarrow). Skipped.")
                return "failed"

            if not self.dry_run:
                meta.setdefault("files", {})
                meta["files"]["cache_umap_obs"] = cache_file
                meta["cache_meta"] = cache_meta
                meta.setdefault("embedding_meta", {})
                meta["embedding_meta"]["cache_umap_obs"] = cache_meta
                meta.setdefault("embedding_files", {})
                meta["embedding_files"]["cache_umap_obs"] = cache_file
                meta.setdefault("embeddings", list(adata.obsm.keys()))

                with open(meta_path, "w") as f:
                    json.dump(meta, f, indent=2)
                print(f"   💾 Cache UMAP+obs: {cache_file}")

            print("   ✅ Done.")
            del adata; gc.collect()
            return "success"
        except Exception as e:
            print(f"   ❌ Cache export failed: {e}")
            return "failed"

    def run_convert(self):
        files = self.scan_r_inputs()
        if not files:
            print("No input files found.")
            return

        print("\n📋 Task List:")
        for i, (s_dir, s_name, _) in enumerate(files, 1):
            print(f"  {i}. {s_dir}/{s_name}")
        print(f"  Total: {len(files)} files\n")

        stats = {'success': 0, 'skipped': 0, 'failed': 0}
        failed_list = []
        total = len(files)

        for i, (s_dir, s_name, fpath) in enumerate(files, 1):
            status = self.process_conversion(s_dir, s_name, fpath, i, total)
            stats[status] += 1
            if status == 'failed':
                failed_list.append(f"{s_dir}/{s_name}")

        print("\n" + "="*60)
        print("🎉 Batch Summary:")
        print(f"  Total   : {len(files)}")
        print(f"  Success : {stats['success']}")
        print(f"  Skipped : {stats['skipped']}")
        print(f"  Failed  : {stats['failed']}")
        
        if failed_list:
            print("\n❌ Failed List:")
            for item in failed_list:
                print(f"  - {item}")
        print("="*60)

    def run_magic(self):
        files = self.scan_h5ad_inputs()
        if not files:
            print("No H5AD files found.")
            return

        print("=" * 60)
        print("✨ MAGIC backfill")
        print(f"Input dir : {self.input_dir}")
        if self.study_filter:
            print(f"Study     : {self.study_filter}")
        print(f"Dry run   : {self.dry_run}")
        print(f"Max cells : {self.magic_params.get('max_cells')}")
        print("=" * 60)

        print("\n📋 Task List:")
        for i, fpath in enumerate(files, 1):
            rel = fpath.relative_to(self.input_dir) if fpath.is_relative_to(self.input_dir) else fpath
            print(f"  {i}. {rel}")
        print(f"  Total: {len(files)} files\n")

        stats = {'success': 0, 'skipped': 0, 'failed': 0}
        total = len(files)
        
        for i, fpath in enumerate(files, 1):
            status = self.process_magic_backfill(fpath, i, total)
            stats[status] += 1

        print("\n" + "="*60)
        print("🎉 Batch Summary:")
        print(f"  Total   : {len(files)}")
        print(f"  Success : {stats['success']}")
        print(f"  Skipped : {stats['skipped']}")
        print(f"  Failed  : {stats['failed']}")
        print("="*60)

    def run_cache(self):
        """Generate cache feather from existing h5ad."""
        files = self.scan_cache_inputs()
        if not files:
            print("No H5AD files found.")
            return

        print("=" * 60)
        print("📦 Cache export (UMAP+obs feather)")
        print(f"Input dir : {self.input_dir}")
        if self.study_filter:
            print(f"Study     : {self.study_filter}")
        print(f"Dry run   : {self.dry_run}")
        print(f"Overwrite : {self.overwrite}")
        print("=" * 60)

        print("\n📋 Task List:")
        for i, fpath in enumerate(files, 1):
            rel = fpath.relative_to(self.input_dir) if fpath.is_relative_to(self.input_dir) else fpath
            print(f"  {i}. {rel}")
        print(f"  Total: {len(files)} files\n")

        stats = {'success': 0, 'skipped': 0, 'failed': 0}
        total = len(files)

        for i, fpath in enumerate(files, 1):
            status = self.process_cache_export(fpath, i, total)
            stats[status] += 1

        print("\n" + "="*60)
        print("🎉 Batch Summary:")
        print(f"  Total   : {len(files)}")
        print(f"  Success : {stats['success']}")
        print(f"  Skipped : {stats['skipped']}")
        print(f"  Failed  : {stats['failed']}")
        print("="*60)

def main():
    parser = argparse.ArgumentParser()
    sub = parser.add_subparsers(dest='cmd', required=True)
    
    # CONVERT
    pc = sub.add_parser('convert', help="Full R -> H5AD conversion")
    pc.add_argument("--input-dir", required=True, type=Path)
    pc.add_argument("--output-dir", required=True, type=Path)
    pc.add_argument("--study", type=str)
    pc.add_argument("--overwrite", action="store_true")
    pc.add_argument("--dry-run", action="store_true")
    pc.add_argument("--force-norm-from", type=str)
    pc.add_argument("--compression", type=str, default="gzip", help="Compression for main h5ad (default gzip; use lzf/none to test speed)")

    # MAGIC
    pm = sub.add_parser('magic', help="Add MAGIC to existing H5AD")
    pm.add_argument("--input-dir", required=True, type=Path)
    pm.add_argument("--study", type=str)
    pm.add_argument("--dry-run", action="store_true")
    pm.add_argument("--max-cells-magic", type=int, default=50000)
    pm.add_argument("--magic-memory-budget", type=float, default=50.0, help="Peak memory budget for MAGIC in GB")

    # CACHE
    pcache = sub.add_parser('cache', help="Export UMAP+obs feather from existing H5AD")
    pcache.add_argument("--input-dir", required=True, type=Path)
    pcache.add_argument("--study", type=str)
    pcache.add_argument("--dry-run", action="store_true")
    pcache.add_argument("--overwrite", action="store_true")

    args = parser.parse_args()
    
    if args.cmd == 'convert':
        ScRNAH5ADConverter(
            input_dir=args.input_dir,
            output_dir=args.output_dir,
            study_filter=args.study,
            overwrite=args.overwrite,
            dry_run=args.dry_run,
            force_norm_from=args.force_norm_from,
            compression=args.compression,
        ).run_convert()
    elif args.cmd == 'magic':
        ScRNAH5ADConverter(
            input_dir=args.input_dir,
            study_filter=args.study,
            magic_params={'max_cells': args.max_cells_magic},
            dry_run=args.dry_run,
            magic_memory_budget=args.magic_memory_budget
        ).run_magic()
    else:
        ScRNAH5ADConverter(
            input_dir=args.input_dir,
            study_filter=args.study,
            dry_run=args.dry_run,
            overwrite=args.overwrite,
        ).run_cache()

if __name__ == "__main__":
    main()
