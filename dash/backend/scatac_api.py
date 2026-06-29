"""
scATAC API (H5AD-based, Gene Activity Score)
Aligns with scRNA API architecture:
- Single h5ad per sample (X = Gene Activity Score, log1p)
- Robust metadata and discovery
- Integrated gene ID resolution (GeneInfo -> CommonID)
"""

from __future__ import annotations

import hashlib
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

logger = logging.getLogger(__name__)


@dataclass(frozen=True)
class _CacheKey:
    dataset_id: str
    mtime: float


class _AnnDataCache:
    """LRU cache for backed AnnData with explicit close."""

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


class scATACAPI:
    def __init__(self, species_key: str):
        if not species_key or not isinstance(species_key, str):
            raise ValueError("Species key is required and must be a non-empty string")
        self.species = species_key.strip()
        if not self.species:
            raise ValueError("Species key cannot be empty or whitespace")

        data_root = os.environ.get("DASH_DATA_PATH", Path(__file__).parent.parent / "data")
        self.data_dir = Path(data_root) / self.species / "scATAC"
        self.logger = logger
        self.logger.info("📂 Loading scATAC studies from %s", self.data_dir)
        
        self.studies = self._discover_studies()
        self._idle_ttl = self._get_int_env("H5AD_IDLE_TTL_SEC", default=900, min_value=60)
        self._idle_sweep = self._get_int_env("H5AD_IDLE_SWEEP_SEC", default=60, min_value=15)
        self._adata_cache = _AnnDataCache(maxsize=16, name="atac")
        self._start_idle_reaper()
        
        # Lazy-loaded helpers
        self._gene_api = None
        self._adapter = None

    # ---------- Discovery & Summary ----------
    def _discover_studies(self) -> List[Dict]:
        studies: List[Dict] = []
        if not self.data_dir.exists():
            self.logger.warning("scATAC data directory not found: %s", self.data_dir)
            return studies

        for study_dir in sorted(self.data_dir.iterdir(), key=lambda p: p.name, reverse=True):
            if not study_dir.is_dir():
                continue
            
            samples = self._get_study_samples(study_dir)
            if samples:
                studies.append({
                    "id": study_dir.name,
                    "name": study_dir.name.replace("_", " "),
                    "path": study_dir,
                    "samples": samples,
                })
        
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

            meta = self._load_metadata(sample_dir, sample_name, h5ad_path)
            
            samples.append({
                "id": f"{study_path.name}/{sample_name}",
                "study_id": study_path.name,
                "sample_name": sample_name,
                "display_name": sample_name.replace("_", " "),
                "path": sample_dir,
                "metadata": meta,
                "usable": True # scATAC validation is simpler for now
            })
        return samples

    def _load_metadata(self, sample_dir: Path, sample_name: str, h5ad_path: Path) -> Dict:
        meta: Dict = {}
        meta_path = sample_dir / f"{sample_name}_metadata.json"
        
        if meta_path.exists():
            try:
                meta = json.loads(meta_path.read_text())
            except Exception as e:
                self.logger.warning(f"Failed to read metadata {meta_path}: {e}")
        
        # Defaults if not in metadata
        meta.setdefault("sample_name", sample_name)
        meta.setdefault("files", {})["h5ad"] = h5ad_path.name
        
        # Check for fragment file
        if "fragment_path" not in meta:
             frag_path = sample_dir / "fragments.tsv.gz"
             if frag_path.exists():
                 meta["fragment_path"] = "fragments.tsv.gz"
        
        return meta

    # ---------- Public API (Listing) ----------
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
                datasets.append({
                    "dataset_id": sample["id"],
                    "study": study["name"],
                    "study_id": study["id"],
                    "sample_name": sample["sample_name"],
                    "name": f"{study['name']} - {sample['display_name']}",
                    "n_cells": meta.get("n_cells", 0),
                    "n_genes": meta.get("n_genes", 0),
                    "cluster_col": meta.get("cluster_col"),
                    "fragment_path": meta.get("fragment_path"),
                    "has_fragment": bool(meta.get("fragment_path")),
                })
        return datasets

    def get_dataset_info(self, dataset_id: str) -> Optional[Dict]:
        for ds in self.get_available_datasets():
            if ds["dataset_id"] == dataset_id:
                return ds
        return None

    # ---------- Data Access (Cached) ----------
    def _resolve_dataset_path(self, dataset_id: str) -> Path:
        if "/" not in dataset_id:
            raise ValueError("dataset_id must be 'study/sample'")
        study, sample = dataset_id.split("/", 1)
        return self.data_dir / study / sample / f"{sample}.h5ad"

    def _get_mtime(self, path: Path) -> float:
        try:
            return path.stat().st_mtime
        except Exception:
            return 0.0

    def _get_adata(self, dataset_id: str) -> sc.AnnData:
        path = self._resolve_dataset_path(dataset_id)
        key = _CacheKey(dataset_id=dataset_id, mtime=self._get_mtime(path))
        removed = self._prune_caches()
        if removed:
            self.logger.info("🧹 h5ad idle clear: %s", ", ".join(removed))
        cached = self._adata_cache.get(key)
        if cached is not None:
            return cached
        if not path.exists():
            raise FileNotFoundError(f"h5ad not found: {path}")
        return self._adata_cache.set(key, lambda: sc.read_h5ad(path, backed="r"))

    def _get_int_env(self, key: str, default: int, min_value: int) -> int:
        try:
            val = int(os.environ.get(key, default))
        except Exception:
            return default
        return max(val, min_value)

    def _prune_caches(self) -> List[str]:
        now = time.time()
        return self._adata_cache.prune_idle(self._idle_ttl, now)

    def _start_idle_reaper(self) -> None:
        def _loop():
            while True:
                time.sleep(self._idle_sweep)
                removed = self._prune_caches()
                if removed:
                    self.logger.info("🧹 h5ad idle clear: %s", ", ".join(removed))

        t = threading.Thread(target=_loop, name="h5ad-idle-reaper", daemon=True)
        t.start()
    
    def _get_sample_entry(self, dataset_id: str) -> Optional[Dict]:
        for study in self.studies:
            for sample in study.get("samples", []):
                if sample.get("id") == dataset_id:
                    return sample
        return None

    def _load_obs_fast(self, sample: Dict) -> Optional[pd.DataFrame]:
        """Prioritize cache/obs feather for metadata, return None on failure."""
        meta = sample.get("metadata", {}) or {}

        # Prioritize cache feather (includes UMAP + obs)
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
                        # Drop embedding columns, keep only obs
                        umap_cols = cache_meta.get("umap_columns") or []
                        df = df.drop(columns=[c for c in umap_cols if c in df.columns], errors="ignore")
                        expected_hash = cache_meta.get("obs_hash")
                        if expected_hash:
                            actual_hash = hashlib.sha1(
                                "\n".join([str(x) for x in df.index]).encode("utf-8")
                            ).hexdigest()
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

        # Fallback to obs feather (for compatibility if needed in future)
        obs_files = meta.get("obs_files") or {}
        if feather is None or not obs_files or "obs" not in obs_files:
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

    # ---------- Metadata & Embeddings ----------
    def get_cell_metadata(self, dataset_id: str) -> pd.DataFrame:
        sample = self._get_sample_entry(dataset_id)
        if sample:
            obs_df = self._load_obs_fast(sample)
            if obs_df is not None:
                return obs_df

        adata = self._get_adata(dataset_id)
        return adata.obs.copy()

    def get_umap_coordinates(self, dataset_id: str) -> Tuple[np.ndarray, np.ndarray]:
        sample = self._get_sample_entry(dataset_id)
        meta = sample.get("metadata", {}) if sample else {}

        # Prioritize cache feather
        cache_meta = meta.get("cache_meta") or {}
        cache_file = (meta.get("files") or {}).get("cache_umap_obs")
        if cache_file and cache_meta and feather is not None and sample:
            fpath = sample["path"] / cache_file
            if fpath.exists():
                try:
                    t0 = time.time()
                    table = feather.read_table(fpath, memory_map=False)
                    if "cell_id" not in table.column_names:
                        raise ValueError(f"cell_id column missing in {fpath}")
                    umap_cols = cache_meta.get("umap_columns") or [
                        c for c in table.column_names if c.startswith("umap_")
                    ]
                    if not umap_cols:
                        raise KeyError("No UMAP columns found in cache feather")
                    coords = np.column_stack([table[c].to_numpy(zero_copy_only=False) for c in umap_cols])
                    cell_ids = table["cell_id"].to_numpy(zero_copy_only=False).astype(str)

                    actual_hash = hashlib.sha1(
                        "\n".join([str(x) for x in cell_ids]).encode("utf-8")
                    ).hexdigest()
                    expected_hash = cache_meta.get("obs_hash")
                    if expected_hash and actual_hash != expected_hash:
                        meta_path = sample.get("path") / (sample.get("sample_name") + "_metadata.json")
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
                    return coords[:, 0], coords[:, 1]
                except Exception as exc:
                    self.logger.warning("Failed to load cache feather %s: %s", fpath, exc)

        adata = self._get_adata(dataset_id)
        
        # Auto-detect embedding
        valid_keys = ["X_umap", "X_tsne", "X_pca", "X_draw_graph_fa"]
        found_key = None
        
        # 1. Check priority keys
        for key in valid_keys:
            if key in adata.obsm_keys():
                found_key = key
                break
        
        # 2. Check any X_ key
        if not found_key:
            for key in adata.obsm_keys():
                if key.startswith("X_") and key != "X_diffmap":
                    found_key = key
                    break
        
        if not found_key:
             raise KeyError(f"No suitable 2D embedding (X_umap, X_tsne, etc.) found in h5ad. Available keys: {list(adata.obsm_keys())}")

        coords = np.asarray(adata.obsm[found_key])
        return coords[:, 0], coords[:, 1]
    
    def get_embeddings(self, dataset_id: str) -> pd.DataFrame:
        x, y = self.get_umap_coordinates(dataset_id)
        return pd.DataFrame({"UMAP1": x, "UMAP2": y}, index=self.get_cell_metadata(dataset_id).index)

    def get_cluster_info(self, dataset_id: str) -> Dict[str, Dict]:
        obs = self.get_cell_metadata(dataset_id)
        if obs.empty:
            return {}
        
        candidates_tier1 = []
        candidates_tier2 = []
        
        # Blocklist
        blocklist = {
            "n_genes", "n_counts", "total_counts", "n_counts_rna", "n_features_rna", 
            "percent_mito", "scrublet_score", "doublet_score", "s_score", "g2m_score",
            "batch", "phase", "barcode", "cell_id"
        }
        
        keywords = ["cluster", "snn", "leiden", "louvain", "celltype", "cell_type", "subtype", "supercelltype"]
        
        # Try to find cluster column from metadata
        sample = self._get_sample_entry(dataset_id)
        meta = sample.get("metadata", {}) if sample else {}
        target_col = meta.get("cluster_col")
        
        # If metadata specifies a column, force include it first
        forced_candidate = None
        if target_col and target_col in obs.columns:
            forced_candidate = target_col
        
        for col in obs.columns:
            if col == forced_candidate:
                continue # Add later manually
                
            col_lower = col.lower()
            if col_lower in blocklist: continue
            if pd.api.types.is_float_dtype(obs[col]): continue
            
            if pd.api.types.is_categorical_dtype(obs[col]):
                num_unique = len(obs[col].cat.categories)
            else:
                num_unique = obs[col].nunique()
            
            if num_unique <= 1 or num_unique > 200: continue
            
            if any(k in col_lower for k in keywords):
                candidates_tier1.append(col)
            else:
                candidates_tier2.append(col)
        
        candidates = candidates_tier1 + candidates_tier2
        if forced_candidate:
            # Metadata specified column goes absolutely first
            if forced_candidate in candidates:
                candidates.remove(forced_candidate)
            candidates.insert(0, forced_candidate)
            
        if not candidates and not obs.columns.empty:
            candidates.append(obs.columns[0])

        cluster_info = {}
        for col in candidates:
             try:
                clusters = obs[col].astype(str)
                counts = clusters.value_counts().sort_index()
                cluster_info[col] = {
                    "n_clusters": len(counts),
                    "cluster_sizes": counts.to_dict(),
                    "total_cells": len(obs),
                }
             except Exception:
                 continue
        return cluster_info

    # ---------- Gene Resolution Helpers ----------
    def _get_gene_api(self):
        if self._gene_api is None:
            try:
                from .gene_search import GeneSearchAPI
                self._gene_api = GeneSearchAPI(self.species)
            except Exception as e:
                self.logger.warning(f"Failed to init GeneSearchAPI: {e}")
                self._gene_api = False
        return self._gene_api if self._gene_api is not False else None

    def _get_adapter(self):
        if self._adapter is None:
            try:
                from .species_adapters.factory import create_species_adapter
                self._adapter = create_species_adapter(self.species)
            except Exception as e:
                self.logger.warning(f"Failed to init adapter: {e}")
                self._adapter = False
        return self._adapter if self._adapter is not False else None

    def _parse_primary_id(self, primary_id: str) -> str:
        if '::' in primary_id:
            return primary_id.split('::', 1)[0]
        return primary_id

    def _resolve_gene_indices(
        self, adata, gene_ids: List[str]
    ) -> Tuple[List[int], Dict[str, str], Dict[str, str]]:
        """
        Resolve primary_ids to var indices for Gene Activity.
        Matches scRNA logic.
        """
        def _norm(val: str) -> Optional[str]:
            if not isinstance(val, str): return None
            v = val.strip()
            return v.lower() if v else None

        # Build lookup from adata.var
        var_names = list(adata.var_names)
        lookup: Dict[str, int] = {}
        for idx, name in enumerate(var_names):
            normed = _norm(name)
            if normed: lookup[normed] = idx
        
        # Check if var has gene_name column
        if "gene_name" in adata.var.columns:
            for idx, val in enumerate(adata.var["gene_name"]):
                normed = _norm(val)
                if normed: lookup[normed] = idx

        # Batch query names from geneinfo
        gene_api = self._get_gene_api()
        id_to_name: Dict[str, str] = {}
        if gene_api:
            try:
                id_to_name = gene_api.get_gene_names_by_ids(gene_ids)
            except Exception:
                pass

        matched_indices: List[int] = []
        resolved_mapping: Dict[str, str] = {}
        display_labels: Dict[str, str] = {}
        seen_hits = set()

        for gid in gene_ids:
            hit = None
            gene_name = id_to_name.get(gid)
            label_source = None

            # 1. Try gene_name
            if gene_name:
                normed = _norm(gene_name)
                if normed and normed in lookup:
                    hit = lookup[normed]
                    label_source = "gene_name"
            
            # 2. Try common_id
            if hit is None:
                common_id = self._parse_primary_id(gid)
                normed = _norm(common_id)
                if normed and normed in lookup:
                    hit = lookup[normed]
                    label_source = "common_id"

            if hit is not None:
                resolved_name = var_names[hit]
                resolved_mapping[gid] = resolved_name
                
                if label_source == "gene_name" and gene_name:
                    display_labels[gid] = resolved_name # Use matrix var name (usually gene name)
                else:
                    common_id = self._parse_primary_id(gid)
                    if gene_name:
                         display_labels[gid] = f"{gene_name} ({common_id})"
                    else:
                         display_labels[gid] = common_id
                
                if hit not in seen_hits:
                    matched_indices.append(hit)
                    seen_hits.add(hit)
        
        return matched_indices, resolved_mapping, display_labels

    # ---------- Gene Activity Access ----------
    def get_gene_activity_batch(
        self, gene_ids: List[str], dataset_id: str
    ) -> Tuple[pd.DataFrame, Dict[str, str], Dict[str, str]]:
        """
        Get Gene Activity Scores for list of gene_ids.
        Returns: (DataFrame [genes x cells], resolved_mapping, display_labels)
        """
        if not gene_ids:
            return pd.DataFrame(), {}, {}
        
        adata = self._get_adata(dataset_id)
        if adata is None:
            raise ValueError(f"Dataset not found: {dataset_id}")
        if adata.X is None:
            raise ValueError("Gene activity matrix (adata.X) is missing in h5ad")
        if adata.n_obs == 0 or adata.n_vars == 0:
            raise ValueError("h5ad has empty cells or genes for Gene Activity")
        matrix = adata.X
        
        matched_indices, resolved_mapping, display_labels = self._resolve_gene_indices(adata, gene_ids)
        
        if not matched_indices:
            return pd.DataFrame(), {}, {}
        
        subset = matrix[:, matched_indices]
        if sp.issparse(subset):
            arr = subset.toarray()
        elif hasattr(subset, "toarray"):
            arr = subset.toarray()
        else:
            arr = np.asarray(subset)
            
        if arr.ndim == 1:
            arr = arr.reshape(-1, 1)
        
        values = arr.T # genes x cells
        
        # Use display labels or var names for index? 
        # scRNA uses var_names from adata. 
        # Let's stick to var_names from adata for consistency in DataFrame, 
        # the caller can use display_labels to rename if needed.
        activity_df = pd.DataFrame(
            values, 
            index=[adata.var_names[i] for i in matched_indices], 
            columns=pd.Index(adata.obs_names, dtype=str)
        )
        
        return activity_df, resolved_mapping, display_labels

    def get_gene_activity_by_cluster(
        self, gene_ids: List[str], cluster_col: str, dataset_id: str
    ) -> pd.DataFrame:
        """
        Aggregate Gene Activity by cluster.
        Returns: DataFrame [genes x clusters]
        """
        activity_df, _, display_labels = self.get_gene_activity_batch(gene_ids, dataset_id)
        if activity_df.empty:
            return pd.DataFrame()
            
        obs = self.get_cell_metadata(dataset_id)
        if obs.empty:
            raise ValueError("Cell metadata is empty.")
        if cluster_col not in obs.columns:
            return pd.DataFrame()
        
        # Align column/index types
        obs_index = pd.Index(obs.index, dtype=str)
        activity_df.columns = pd.Index(activity_df.columns, dtype=str)

        common_cells = activity_df.columns.intersection(obs_index)
        if common_cells.empty:
            return pd.DataFrame()
            
        subset = activity_df[common_cells]
        clusters = obs.loc[common_cells, cluster_col].astype(str)
        
        subset_T = subset.T.copy()
        subset_T["__cluster__"] = clusters.values
        grouped = subset_T.groupby("__cluster__").mean().T
        grouped.index.name = "Gene"
        # Prefer user-facing labels when available
        if display_labels:
            grouped = grouped.rename(index=display_labels)
        return grouped
