"""
H5AD-based scRNA API (replaces legacy multi-file implementation).

Design:
- Single h5ad file: X=norm, layers may contain raw/magic, obsm contains embeddings.
- Provides same main interfaces as legacy scRNAAPI for seamless Dash frontend usage.
"""

from __future__ import annotations

import json
import logging
import os
import threading
import time
from dataclasses import dataclass
from collections import OrderedDict
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp
try:
    import pyarrow.feather as feather
except ImportError:
    feather = None
import hashlib

logger = logging.getLogger(__name__)


@dataclass(frozen=True)
class _CacheKey:
    dataset_id: str
    mtime: float


class _AnnDataCache:
    """LRU cache for backed AnnData objects with explicit close on eviction."""

    def __init__(self, maxsize: int, name: str):
        self.maxsize = maxsize
        self.name = name
        self._store: "OrderedDict[_CacheKey, Tuple[Any, float]]" = OrderedDict()
        self._lock = threading.Lock()

    def _now(self) -> float:
        return time.time()

    def _close(self, adata: Any) -> None:
        try:
            file_handle = getattr(adata, "file", None)
            if file_handle is not None:
                file_handle.close()
        except Exception:
            pass

    def get(self, key: _CacheKey):
        with self._lock:
            entry = self._store.get(key)
            if entry is None:
                return None
            adata, _ = entry
            self._store.move_to_end(key)
            self._store[key] = (adata, self._now())
            return adata

    def set(self, key: _CacheKey, loader) -> Any:
        with self._lock:
            entry = self._store.get(key)
            if entry is not None:
                adata, _ = entry
                self._store.move_to_end(key)
                self._store[key] = (adata, self._now())
                return adata
        adata = loader()
        with self._lock:
            self._store[key] = (adata, self._now())
            self._evict_if_needed()
        return adata

    def _evict_if_needed(self) -> None:
        while len(self._store) > self.maxsize:
            _, (adata, _) = self._store.popitem(last=False)
            self._close(adata)

    def prune_idle(self, ttl_seconds: int, now: Optional[float] = None) -> List[str]:
        if ttl_seconds <= 0:
            return []
        now_ts = now or self._now()
        removed: List[str] = []
        with self._lock:
            for key, (adata, last_use) in list(self._store.items()):
                if now_ts - last_use < ttl_seconds:
                    break
                self._store.pop(key, None)
                removed.append(f"{self.name}:{key.dataset_id}")
                self._close(adata)
        return removed

    def clear(self) -> None:
        with self._lock:
            for adata, _ in self._store.values():
                self._close(adata)
            self._store.clear()


class scRNAAPI:
    def __init__(self, species_key: str):
        if not species_key or not isinstance(species_key, str):
            raise ValueError("Species key is required and must be a non-empty string")
        self.species = species_key.strip()
        if not self.species:
            raise ValueError("Species key cannot be empty or whitespace")

        data_root = os.environ.get("DASH_DATA_PATH", Path(__file__).parent.parent / "data")
        self.data_dir = Path(data_root) / self.species / "scRNA"
        self.logger = logger
        self.logger.info("📂 Loading scRNA studies from %s", self.data_dir)
        self.studies = self._discover_studies()
        # Cache usability map for quick guard checks
        self._dataset_usable = {
            sample["id"]: sample.get("usable", True)
            for study in self.studies
            for sample in study.get("samples", [])
        }
        self._idle_ttl = self._get_int_env("H5AD_IDLE_TTL_SEC", default=900, min_value=60)
        self._idle_sweep = self._get_int_env("H5AD_IDLE_SWEEP_SEC", default=60, min_value=15)
        self._adata_cache = _AnnDataCache(maxsize=32, name="adata")
        self._magic_cache = _AnnDataCache(maxsize=16, name="magic")
        self._start_idle_reaper()

        # Lazy-loaded helpers for gene ID resolution
        self._gene_api = None
        self._adapter = None

    # ---------- discovery & summary ----------
    def _discover_studies(self) -> List[Dict]:
        studies: List[Dict] = []
        if not self.data_dir.exists():
            self.logger.warning("scRNA data directory not found: %s", self.data_dir)
            return studies

        total_samples = 0
        invalid_samples = 0
        total_start = time.time()
        for study_dir in sorted(self.data_dir.iterdir(), key=lambda p: p.name, reverse=True):
            if not study_dir.is_dir():
                continue
            study_start = time.time()
            samples = self._get_study_samples(study_dir)
            if not samples:
                continue
            total_samples += len(samples)
            invalid_samples += sum(1 for s in samples if not s.get("usable", True))
            self.logger.info(
                "📖 Study %s: %d samples (invalid %d) in %.2fs",
                study_dir.name,
                len(samples),
                sum(1 for s in samples if not s.get("usable", True)),
                time.time() - study_start,
            )
            studies.append(
                {
                    "id": study_dir.name,
                    "name": study_dir.name.replace("_", " "),
                    "path": study_dir,
                    "samples": samples,
                }
            )
        self.logger.info(
            "✅ scRNA discovery done: %d studies, %d samples (invalid %d) in %.2fs",
            len(studies),
            total_samples,
            invalid_samples,
            time.time() - total_start,
        )
        return studies

    def _get_study_samples(self, study_path: Path) -> List[Dict]:
        samples: List[Dict] = []
        for sample_dir in sorted(study_path.iterdir()):
            if not sample_dir.is_dir():
                continue
            sample_name = sample_dir.name
            h5ad_path = sample_dir / f"{sample_name}.h5ad"
            if not h5ad_path.exists():
                continue
            t0 = time.time()
            meta = self._load_manifest_or_summary(h5ad_path)
            usable, reason = self._validate_metadata(meta)
            self.logger.info(
                "• %s/%s -> %s in %.2fs%s",
                study_path.name,
                sample_name,
                "usable" if usable else "UNUSABLE",
                time.time() - t0,
                f" ({reason})" if not usable else "",
            )
            samples.append(
                {
                    "id": f"{study_path.name}/{sample_name}",
                    "study_id": study_path.name,
                    "sample_name": sample_name,
                    "display_name": sample_name.replace("_", " "),
                    "path": sample_dir,
                    "metadata": meta,
                    "usable": usable,
                    "unusable_reason": reason,
                }
            )
        return samples

    def _load_manifest_or_summary(self, h5ad_path: Path) -> Dict:
        meta: Dict = {}
        manifest = h5ad_path.with_name(f"{h5ad_path.stem}_metadata.json")
        if manifest.exists():
            try:
                meta = json.loads(manifest.read_text())
            except Exception as exc:
                self.logger.warning("Failed to read manifest %s: %s", manifest, exc)

        # If manifest is missing, probe h5ad once to populate summary.
        if not meta:
            try:
                adata = sc.read_h5ad(h5ad_path, backed="r")
                meta = {
                    "n_cells": adata.n_obs,
                    "n_genes": adata.n_vars,
                    "available_data_versions": self._available_versions_from_uns(adata),
                    "default_data_version": self._default_version_from_uns(adata),
                    "embeddings": list(adata.obsm.keys()),
                    "cell_metadata_columns": list(adata.obs.columns),
                    "files": {"h5ad": h5ad_path.name},
                }
                adata.file.close()
            except Exception as exc:
                self.logger.error("Error reading summary from %s: %s", h5ad_path, exc)
        meta.setdefault("assay_used", "RNA")
        meta.setdefault("available_assays", ["RNA"])
        return meta

    def _validate_metadata(self, meta: Dict) -> Tuple[bool, Optional[str]]:
        required_keys = ["n_cells", "n_genes", "available_data_versions", "default_data_version", "files"]
        missing = [k for k in required_keys if k not in meta]
        if missing:
            return False, f"missing fields: {missing}"
        if not meta.get("available_data_versions"):
            return False, "available_data_versions empty"
        if not meta.get("default_data_version"):
            return False, "default_data_version empty"
        if meta.get("n_cells", 0) <= 0:
            return False, "n_cells <= 0"
        if meta.get("n_genes", 0) <= 0:
            return False, "n_genes <= 0"
        files = meta.get("files") or {}
        if "h5ad" not in files:
            return False, "files.h5ad missing"
        return True, None

    def _resolve_dataset_path(self, dataset_id: str, magic: bool = False) -> Path:
        if "/" not in dataset_id:
            raise ValueError("dataset_id must be 'study/sample'")
        study, sample = dataset_id.split("/", 1)
        fname = f"{sample}_magic.h5ad" if magic else f"{sample}.h5ad"
        return self.data_dir / study / sample / fname

    def _get_sample_entry(self, dataset_id: str) -> Optional[Dict]:
        for study in self.studies:
            for sample in study.get("samples", []):
                if sample.get("id") == dataset_id:
                    return sample
        return None

    def _get_mtime(self, path: Path) -> float:
        try:
            return path.stat().st_mtime
        except Exception:
            return 0.0

    def _get_adata(self, dataset_id: str):
        removed = self._prune_caches()
        if removed:
            self.logger.info("🧹 h5ad idle clear: %s", ", ".join(removed))
        if not self._dataset_usable.get(dataset_id, True):
            raise ValueError(f"Dataset {dataset_id} is marked unusable by metadata validation.")
        path = self._resolve_dataset_path(dataset_id)
        key = _CacheKey(dataset_id=dataset_id, mtime=self._get_mtime(path))
        cached = self._adata_cache.get(key)
        if cached is not None:
            return cached
        if not path.exists():
            raise FileNotFoundError(f"h5ad not found: {path}")
        return self._adata_cache.set(key, lambda: sc.read_h5ad(path, backed="r"))

    def _get_magic_adata(self, dataset_id: str):
        removed = self._prune_caches()
        if removed:
            self.logger.info("🧹 h5ad idle clear: %s", ", ".join(removed))
        path = self._resolve_dataset_path(dataset_id, magic=True)
        key = _CacheKey(dataset_id=dataset_id, mtime=self._get_mtime(path))
        cached = self._magic_cache.get(key)
        if cached is not None:
            return cached
        if not path.exists():
            raise FileNotFoundError(f"magic h5ad not found: {path}")
        return self._magic_cache.set(key, lambda: sc.read_h5ad(path, backed="r"))

    def _get_int_env(self, key: str, default: int, min_value: int) -> int:
        try:
            val = int(os.environ.get(key, default))
        except Exception:
            return default
        return max(val, min_value)

    def _prune_caches(self) -> List[str]:
        now = time.time()
        removed = []
        removed.extend(self._adata_cache.prune_idle(self._idle_ttl, now))
        removed.extend(self._magic_cache.prune_idle(self._idle_ttl, now))
        return removed

    def _start_idle_reaper(self) -> None:
        def _loop():
            while True:
                time.sleep(self._idle_sweep)
                removed = self._prune_caches()
                if removed:
                    self.logger.info("🧹 h5ad idle clear: %s", ", ".join(removed))

        t = threading.Thread(target=_loop, name="h5ad-idle-reaper", daemon=True)
        t.start()

    def _load_obs_fast(self, sample: Dict) -> Optional[pd.DataFrame]:
        """Load obs from feather/combined if available; return None on failure."""
        meta = sample.get("metadata", {}) or {}

        # Prefer cache feather (umap+obs) if present
        cache_meta = meta.get("cache_meta") or {}
        cache_file = (meta.get("files") or {}).get("cache_umap_obs")
        if cache_file and cache_meta and feather is not None:
            fpath = sample["path"] / cache_file
            if fpath.exists():
                try:
                    t0 = time.time()
                    table = feather.read_table(fpath, memory_map=False)
                    if "cell_id" in table.column_names:
                        df = table.to_pandas(types_mapper=None).set_index("cell_id")
                        # drop umap columns if present
                        umap_cols = cache_meta.get("umap_columns") or []
                        df = df.drop(columns=[c for c in umap_cols if c in df.columns], errors="ignore")
                        expected_hash = cache_meta.get("obs_hash")
                        if expected_hash:
                            actual_hash = hashlib.sha1("\n".join([str(x) for x in df.index]).encode("utf-8")).hexdigest()
                            if actual_hash != expected_hash:
                                self.logger.warning(
                                    "cache feather hash mismatch (got %s, expected %s) for %s; fallback to h5ad",
                                    actual_hash,
                                    expected_hash,
                                    fpath,
                                )
                            else:
                                self.logger.info(
                                    "⏱️ Loaded cache feather(obs) %s in %.2fs (%s rows, %s cols)",
                                    fpath.name,
                                    time.time() - t0,
                                    df.shape[0],
                                    df.shape[1],
                                )
                                return df
                except Exception as exc:
                    self.logger.warning("Failed to load cache feather %s: %s", fpath, exc)

        obs_files = meta.get("obs_files") or {}
        if not obs_files or "obs" not in obs_files:
            return None
        if feather is None:
            return None
        fpath = sample["path"] / obs_files["obs"]
        if not fpath.exists():
            return None
        try:
            t0 = time.time()
            table = feather.read_table(fpath, memory_map=False)
            if "cell_id" not in table.column_names:
                return None
            df = table.to_pandas(types_mapper=None).set_index("cell_id")
            expected_hash = None
            obs_meta = meta.get("obs_meta") or {}
            if isinstance(obs_meta, dict):
                expected_hash = obs_meta.get("obs", {}).get("obs_hash")
            if expected_hash:
                actual_hash = hashlib.sha1("\n".join([str(x) for x in df.index]).encode("utf-8")).hexdigest()
                if actual_hash != expected_hash:
                    self.logger.warning(
                        "obs feather hash mismatch (got %s, expected %s) for %s; fallback to h5ad",
                        actual_hash,
                        expected_hash,
                        fpath,
                    )
                    return None
            self.logger.info(
                "⏱️ Loaded obs feather %s in %.2fs (%s rows, %s cols)",
                fpath.name,
                time.time() - t0,
                df.shape[0],
                df.shape[1],
            )
            return df
        except Exception as exc:
            self.logger.warning("Failed to load obs feather %s: %s", fpath, exc)
            return None

    # ---------- public API (listing) ----------
    def get_available_studies(self) -> List[Dict]:
        return self.studies

    def get_samples_by_study(self, study_id: str) -> List[Dict]:
        for study in self.studies:
            if study["id"] == study_id:
                return study["samples"]
        return []

    def get_available_datasets(self) -> List[Dict]:
        datasets = []
        for study in self.studies:
            for sample in study["samples"]:
                meta = sample.get("metadata", {})
                available_versions = meta.get("available_data_versions", {})
                default_version = meta.get("default_data_version", "norm")
                usable = sample.get("usable", True)
                datasets.append(
                    {
                        "dataset_id": sample["id"],
                        "study": study["name"],
                        "study_id": study["id"],
                        "sample_name": sample["sample_name"],
                        "name": f"{study['name']} - {sample['display_name']}",
                        "description": meta.get("description", f"scRNA-seq data for {sample['display_name']}"),
                        "n_cells": meta.get("n_cells", 0),
                        "n_genes": meta.get("n_genes", 0),
                        "embeddings": meta.get("embeddings", []),
                        "cell_metadata_columns": meta.get("cell_metadata_columns", []),
                        "available_data_versions": available_versions,
                        "default_data_version": default_version,
                        "has_magic": "magic" in available_versions,
                        "assay_used": meta.get("assay_used", "RNA"),
                        "available_assays": meta.get("available_assays", ["RNA"]),
                        "data_type": meta.get("data_type", "unknown"),
                        "usable": usable,
                        "unusable_reason": sample.get("unusable_reason"),
                    }
                )
        return datasets

    def get_dataset_info(self, dataset_id: str) -> Optional[Dict]:
        for ds in self.get_available_datasets():
            if ds["dataset_id"] == dataset_id:
                return ds
        return None

    # ---------- version selection ----------
    def get_available_data_versions(self, dataset_id: str) -> Dict[str, Dict]:
        sample = self._get_sample_entry(dataset_id)
        if sample:
            meta = sample.get("metadata", {}) or {}
            avail = meta.get("available_data_versions")
            if isinstance(avail, dict) and avail:
                return avail
        # fallback to reading h5ad
        adata = self._get_adata(dataset_id)
        return self._available_versions_from_uns(adata)

    def determine_data_version(
        self, dataset_id: str, requested_version: str = "norm", available_versions: Optional[Dict[str, Dict]] = None
    ) -> Tuple[str, Dict[str, str]]:
        if available_versions is None:
            available_versions = self.get_available_data_versions(dataset_id) or {}
        if not available_versions:
            # last resort: read from h5ad
            adata = self._get_adata(dataset_id)
            available_versions = self._available_versions_from_uns(adata) or {}

        selection = (requested_version or "norm").lower()
        if selection == "auto":
            # try metadata default first
            sample = self._get_sample_entry(dataset_id)
            if sample:
                meta = sample.get("metadata", {}) or {}
                selection = meta.get("default_data_version", "norm")
            else:
                selection = self._default_version_from_uns(self._get_adata(dataset_id))
        if selection not in {"norm", "magic", "raw"}:
            selection = "norm"
        order = {
            "magic": ["magic", "norm", "raw"],
            "norm": ["norm", "magic", "raw"],
            "raw": ["raw", "norm", "magic"],
        }
        resolved = None
        for cand in order[selection]:
            if cand in available_versions:
                resolved = cand
                break
        if resolved is None:
            resolved = selection
        return resolved, available_versions.get(resolved, {})

    # ---------- embeddings & metadata ----------
    def get_embedding_coordinates(self, dataset_id: str) -> Tuple[np.ndarray, np.ndarray, str]:
        sample = self._get_sample_entry(dataset_id)
        if not sample:
            raise KeyError(f"dataset not found: {dataset_id}")
        meta = sample.get("metadata", {}) or {}

        # Prefer cache feather if present
        cache_meta = meta.get("cache_meta") or {}
        cache_file = (meta.get("files") or {}).get("cache_umap_obs")
        if cache_file and cache_meta and feather is not None:
            fpath = sample["path"] / cache_file
            if fpath.exists():
                t0 = time.time()
                table = feather.read_table(fpath, memory_map=False)
                if "cell_id" not in table.column_names:
                    raise ValueError(f"cell_id column missing in {fpath}")
                umap_cols = cache_meta.get("umap_columns") or [c for c in table.column_names if c.startswith("umap_")]
                if not umap_cols:
                    raise KeyError("No UMAP columns found in cache feather")
                coords = np.column_stack([table[c].to_numpy(zero_copy_only=False) for c in umap_cols])
                cell_ids = table["cell_id"].to_numpy(zero_copy_only=False).astype(str)

                actual_hash = hashlib.sha1("\n".join([str(x) for x in cell_ids]).encode("utf-8")).hexdigest()
                expected_hash = cache_meta.get("obs_hash")
                if expected_hash and actual_hash != expected_hash:
                    meta_path = sample.get('path') / (sample.get('sample_name') + '_metadata.json')
                    raise ValueError(
                        "Embedding cell_id hash mismatch metadata.\n"
                        f"  dataset : {dataset_id}\n"
                        f"  got     : {actual_hash}\n"
                        f"  expected: {expected_hash}\n"
                        f"  metadata: {meta_path}\n"
                        f"  embedding: {fpath}\n"
                        "  hint    : rebuild embedding cache"
                    )
                self.logger.info(
                    "⏱️ Loaded cache feather %s in %.2fs (%s rows)",
                    fpath.name,
                    time.time() - t0,
                    coords.shape[0],
                )
                return coords[:, 0], coords[:, 1], "UMAP"

        # Final fallback: read from h5ad obsm
        adata = self._get_adata(dataset_id)
        
        # Auto-detect embedding
        valid_keys = ["X_umap", "X_tsne", "X_pca", "X_draw_graph_fa"]
        found_key = None
        method_name = "UMAP"
        
        # 1. Check priority keys
        for key in valid_keys:
            if key in adata.obsm_keys():
                found_key = key
                method_name = key.replace("X_", "").replace("draw_graph_fa", "ForceAtlas").upper()
                if method_name == "TSNE": method_name = "tSNE"
                break
        
        # 2. Check any X_ key
        if not found_key:
            for key in adata.obsm_keys():
                if key.startswith("X_") and key != "X_diffmap":  # diffmap usually not for 2D vis
                    found_key = key
                    method_name = key.replace("X_", "").upper()
                    break
        
        if not found_key:
            raise KeyError(f"No suitable 2D embedding (X_umap, X_tsne, etc.) found in h5ad. Available keys: {list(adata.obsm_keys())}")

        arr = np.array(adata.obsm[found_key])
        if arr.shape[1] < 2:
             raise ValueError(f"Embedding {found_key} has fewer than 2 dimensions")

        coords = arr[:, :2]
        self.logger.info(
            "⏱️ Loaded %s (%s) from h5ad obsm (no feather) (%s rows)",
            method_name,
            found_key,
            coords.shape[0],
        )
        return coords[:, 0], coords[:, 1], method_name

    def get_umap_coordinates(self, dataset_id: str) -> Tuple[np.ndarray, np.ndarray]:
        x, y, _ = self.get_embedding_coordinates(dataset_id)
        return x, y

    def get_embeddings(self, dataset_id: str) -> pd.DataFrame:
        adata = self._get_adata(dataset_id)
        # Flatten obsm into a dataframe with column names like "{key}_{i}"
        frames = []
        for key in adata.obsm.keys():
            arr = np.array(adata.obsm[key])
            if arr.ndim != 2:
                continue
            cols = [f"{key}_{i+1}" for i in range(arr.shape[1])]
            frames.append(pd.DataFrame(arr, index=adata.obs_names, columns=cols))
        if not frames:
            return pd.DataFrame(index=adata.obs_names)
        return pd.concat(frames, axis=1)

    def get_cell_metadata(self, dataset_id: str) -> pd.DataFrame:
        sample = self._get_sample_entry(dataset_id)
        if not sample:
            raise KeyError(f"dataset not found: {dataset_id}")

        # Prioritize feather obs
        obs_df = self._load_obs_fast(sample)
        if obs_df is not None:
            return obs_df

        adata = self._get_adata(dataset_id)
        return adata.obs.copy()

    def get_gene_metadata(self, dataset_id: str) -> pd.DataFrame:
        adata = self._get_adata(dataset_id)
        return adata.var.copy()

    def get_cluster_info(self, dataset_id: str) -> Dict[str, any]:
        obs = self.get_cell_metadata(dataset_id)
        if obs.empty:
            return {}

        candidates_tier1 = []  # Contains keywords
        candidates_tier2 = []  # Valid but no keywords
        
        # Blacklist (lowercase) - common QC metrics and IDs
        blocklist = {
            "n_genes", "n_counts", "total_counts", "n_counts_rna", "n_features_rna", 
            "percent_mito", "scrublet_score", "doublet_score", "s_score", "g2m_score",
            "batch", "phase", "barcode", "cell_id"
        }
        
        # Priority keywords
        keywords = ["cluster", "snn", "leiden", "louvain", "celltype", "cell_type", "subtype", "supercelltype"]

        for col in obs.columns:
            col_lower = col.lower()
            
            # 1. Blacklist check
            if col_lower in blocklist:
                continue
            
            # 2. Type check: exclude float
            if pd.api.types.is_float_dtype(obs[col]):
                continue

            # 3. Cardinality check
            if pd.api.types.is_categorical_dtype(obs[col]):
                num_unique = len(obs[col].cat.categories)
            else:
                num_unique = obs[col].nunique()
            
            # Rule A: Must have more than 1 group
            if num_unique <= 1:
                continue
                
            # Rule B: Too many groups imply ID or continuous-like data
            if num_unique > 200: 
                continue

            # 4. Classification
            if any(k in col_lower for k in keywords):
                candidates_tier1.append(col)
            else:
                candidates_tier2.append(col)

        # Merge tiers
        candidates = candidates_tier1 + candidates_tier2
        
        # Fallback: if absolutely nothing found, try first column (legacy behavior)
        if not candidates and len(obs.columns) > 0:
            candidates = [obs.columns[0]]

        cluster_info: Dict[str, Dict] = {}
        for col in candidates:
            try:
                clusters = obs[col].astype(str)
                cluster_counts = clusters.value_counts().sort_index()
                cluster_info[col] = {
                    "n_clusters": len(cluster_counts),
                    "cluster_sizes": cluster_counts.to_dict(),
                    "total_cells": len(clusters),
                }
            except Exception as e:
                self.logger.warning("Failed to process cluster column %s: %s", col, e)
        return cluster_info

    # ---------- gene resolution helpers ----------
    def _get_gene_api(self):
        """Lazy-load gene_search_api for ID→name resolution"""
        if self._gene_api is None:
            try:
                from .gene_search import GeneSearchAPI
                self._gene_api = GeneSearchAPI(self.species)
            except Exception as e:
                self.logger.warning(f"Failed to init GeneSearchAPI: {e}")
                self._gene_api = False  # Mark as failed to avoid retry
        return self._gene_api if self._gene_api is not False else None

    def _get_adapter(self):
        """Lazy-load species adapter for parse_primary_id"""
        if self._adapter is None:
            try:
                from .species_adapters.factory import create_species_adapter
                self._adapter = create_species_adapter(self.species)
            except Exception as e:
                self.logger.warning(f"Failed to init adapter: {e}")
                self._adapter = False
        return self._adapter if self._adapter is not False else None

    def _parse_primary_id(self, primary_id: str) -> str:
        """Extract common_id from primary_id (e.g., 'IGDB::ENS' -> 'IGDB')"""
        if '::' in primary_id:
            return primary_id.split('::', 1)[0]
        return primary_id

    # ---------- gene resolution ----------
    def _available_versions_from_uns(self, adata) -> Dict:
        uns = getattr(adata, "uns", {}) or {}
        avail = {}
        if "available_versions" in uns:
            avail = dict(uns["available_versions"])
        else:
            avail["norm"] = {"type": "normalized", "description": "log1p normalized"}
            if "raw" in adata.layers:
                avail["raw"] = {"type": "raw_counts", "description": "raw counts"}
            if "magic" in adata.layers:
                avail["magic"] = {"type": "magic_imputed", "description": "MAGIC smoothed"}
        return avail

    def _default_version_from_uns(self, adata) -> str:
        uns = getattr(adata, "uns", {}) or {}
        if "default_version" in uns:
            return str(uns["default_version"])
        return "norm"

    def _get_matrix(self, adata, layer: str):
        if layer in {"norm", "x", "X"}:
            return adata.X
        if layer in adata.layers:
            return adata.layers[layer]
        return adata.X

    def _resolve_gene_indices(
        self, adata, gene_ids: List[str]
    ) -> Tuple[List[int], Dict[str, str], Dict[str, str]]:
        """
        Resolve primary_ids to var indices (following gene_id_system.md):
        1) primary_id -> gene_name (from geneinfo)
        2) primary_id -> common_id (parse_primary_id, handles IGDB::ENS)
        No "direct input" fallback to ensure consistent primary_id resolution.
        """
        def _norm(val: str) -> Optional[str]:
            if not isinstance(val, str):
                return None
            v = val.strip()
            if not v:
                return None
            return v.lower()

        # Build h5ad index: var_names + var columns (gene_id, gene_name if present)
        var_names = list(adata.var_names)
        lookup: Dict[str, int] = {}
        for idx, name in enumerate(var_names):
            normed = _norm(name)
            if normed:
                lookup[normed] = idx

        for col_name in ("gene_id", "gene_name"):
            if col_name in adata.var.columns:
                for idx, val in enumerate(adata.var[col_name]):
                    normed = _norm(val)
                    if normed:
                        lookup[normed] = idx

        # Batch query: primary_id → gene_name (from geneinfo.db)
        gene_api = self._get_gene_api()
        id_to_name: Dict[str, str] = {}
        if gene_api:
            try:
                id_to_name = gene_api.get_gene_names_by_ids(gene_ids)
            except Exception as e:
                self.logger.warning(f"Failed to lookup gene names: {e}")

        # Dual-lookup strategy (gene_name first, then common_id)
        matched_indices: List[int] = []
        resolved_mapping: Dict[str, str] = {}
        display_labels: Dict[str, str] = {}
        seen_hits = set()

        for gid in gene_ids:
            hit = None

            gene_name = id_to_name.get(gid)
            label_source = None

            # Strategy 1: gene_name
            if gene_name:
                normed = _norm(gene_name)
                if normed and normed in lookup:
                    hit = lookup[normed]
                    label_source = "gene_name"

            # Strategy 2: common_id
            if hit is None:
                common_id = self._parse_primary_id(gid)
                normed = _norm(common_id)
                if normed and normed in lookup:
                    hit = lookup[normed]
                    label_source = "common_id"

            # Record match
            if hit is not None:
                resolved_name = var_names[hit]
                resolved_mapping[gid] = resolved_name

                if label_source == "gene_name" and gene_name:
                    # Use dataset's row name to avoid case confusion
                    display_labels[gid] = resolved_name
                else:
                    common_id = self._parse_primary_id(gid)
                    # Matched by ID: show "gene_name (common_id)" if gene_name exists, else common_id
                    if gene_name:
                        display_labels[gid] = f"{gene_name} ({common_id})"
                    else:
                        display_labels[gid] = common_id

                if hit not in seen_hits:
                    matched_indices.append(hit)
                    seen_hits.add(hit)
            else:
                self.logger.debug(
                    "Gene %s not found in dataset (tried gene_name=%s, common_id=%s)",
                    gid,
                    gene_name,
                    self._parse_primary_id(gid),
                )

        return matched_indices, resolved_mapping, display_labels

    # ---------- expression ----------
    def get_gene_expression_batch(
        self, gene_ids: List[str], dataset_id: str, data_version: str = "norm"
    ) -> Tuple[pd.DataFrame, Dict[str, str], Dict[str, str]]:
        if not gene_ids:
            return pd.DataFrame(), {}, {}
        sample = self._get_sample_entry(dataset_id)
        meta = sample.get("metadata", {}) if sample else {}

        # Decide which file to read based on resolved layer (prefer metadata to avoid hitting h5ad)
        available_versions = meta.get("available_data_versions") or {}
        resolved_layer, _ = self.determine_data_version(dataset_id, data_version, available_versions)
        adata_main = None

        use_magic = resolved_layer == "magic"
        if use_magic:
            magic_file = (meta.get("files") or {}).get("magic_h5ad")
            if not magic_file:
                raise ValueError("magic requested but magic_h5ad not available")
            adata = self._get_magic_adata(dataset_id)
        else:
            if adata_main is None:
                adata_main = self._get_adata(dataset_id)
            adata = adata_main
        matrix = self._get_matrix(adata, resolved_layer)

        matched_indices, resolved_mapping, display_labels = self._resolve_gene_indices(adata, gene_ids)
        if not matched_indices:
            return pd.DataFrame(), {}, {}

        if matrix.shape[1] < max(matched_indices) + 1:
            raise IndexError("Gene index out of bounds for matrix")

        subset = matrix[:, matched_indices]
        if sp.issparse(subset):
            arr = subset.toarray()
        elif hasattr(subset, "toarray"):
            arr = subset.toarray()
        else:
            arr = np.asarray(subset)
        if arr.ndim == 1:
            arr = arr.reshape(-1, 1)
        values = arr.T  # genes × cells

        expr_df = pd.DataFrame(values, index=[adata.var_names[i] for i in matched_indices], columns=adata.obs_names)
        return expr_df, resolved_mapping, display_labels

    def get_gene_expression(self, gene_ids: List[str], dataset_id: str, data_version: str = "norm") -> pd.DataFrame:
        expr_df, _, _ = self.get_gene_expression_batch(gene_ids, dataset_id, data_version)
        return expr_df

    def get_expression_by_cluster(
        self,
        gene_ids: List[str],
        cluster_column: str,
        dataset_id: str,
        data_version: str = "norm",
    ) -> pd.DataFrame:
        expr_df, _, _ = self.get_gene_expression_batch(gene_ids, dataset_id, data_version)
        if expr_df.empty:
            return pd.DataFrame()
        obs = self.get_cell_metadata(dataset_id)
        if obs.empty or cluster_column not in obs.columns:
            return pd.DataFrame()

        common_cells = expr_df.columns.intersection(obs.index)
        if common_cells.empty:
            return pd.DataFrame()
        expr_subset = expr_df[common_cells]
        clusters = obs.loc[common_cells, cluster_column].astype(str)

        expr_with_cluster = expr_subset.T.copy()
        expr_with_cluster["__cluster__"] = clusters.values
        grouped = expr_with_cluster.groupby("__cluster__").mean().T
        grouped.index.name = "Gene"
        grouped = grouped.rename(columns=lambda c: f"Cluster_{c}")
        return grouped
