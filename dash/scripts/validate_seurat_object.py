#!/usr/bin/env python3
"""
Validate Seurat objects against scRNA or scATAC conversion requirements.

This script checks whether a Seurat object (.rds/.Rdata) meets the requirements
specified in dash/docs/data_managers/scrna_conversion.md and scatac_conversion.md.
"""

import argparse
import sys
import os
from pathlib import Path
from typing import Dict, List, Tuple, Optional
import json

try:
    import rpy2.robjects as ro
    from rpy2.robjects import pandas2ri, numpy2ri
    from rpy2.robjects.conversion import localconverter
    from rpy2.rinterface_lib import callbacks
except ImportError:
    print("Error: rpy2 is required. Install with: pip install rpy2", file=sys.stderr)
    sys.exit(1)

try:
    import numpy as np
except ImportError:
    print("Error: numpy is required. Install with: pip install numpy", file=sys.stderr)
    sys.exit(1)


# Suppress R output by default
_suppress_r_output = True


def _custom_write_console(s):
    """Custom callback to suppress R console output"""
    if _suppress_r_output:
        # Suppress all R output including reticulate warnings
        return
    sys.stderr.write(s)


def setup_r_environment(quiet: bool = True):
    """Setup R environment and suppress startup messages"""
    global _suppress_r_output
    _suppress_r_output = quiet

    if quiet:
        # Set custom write console callback to suppress output
        callbacks.consolewrite_print = _custom_write_console
        callbacks.consolewrite_warnerror = _custom_write_console



class ValidationResult:
    """Store validation results"""
    def __init__(self):
        self.passed: List[str] = []
        self.warnings: List[str] = []
        self.errors: List[str] = []
        self.info: Dict[str, any] = {}

    def add_pass(self, message: str):
        self.passed.append(message)

    def add_warning(self, message: str):
        self.warnings.append(message)

    def add_error(self, message: str):
        self.errors.append(message)

    def is_valid(self) -> bool:
        return len(self.errors) == 0

    def print_summary(self, verbose: bool = False):
        """Print validation summary"""
        print("\n" + "="*60)
        print("VALIDATION SUMMARY")
        print("="*60)

        if self.info:
            print("\nBasic Information:")
            for key, value in self.info.items():
                print(f"  {key}: {value}")

        if verbose and self.passed:
            print(f"\n✓ Passed ({len(self.passed)}):")
            for msg in self.passed:
                print(f"  ✓ {msg}")

        if self.warnings:
            print(f"\n⚠ Warnings ({len(self.warnings)}):")
            for msg in self.warnings:
                print(f"  ⚠ {msg}")

        if self.errors:
            print(f"\n✗ Errors ({len(self.errors)}):")
            for msg in self.errors:
                print(f"  ✗ {msg}")

        print("\n" + "="*60)
        if self.is_valid():
            print("RESULT: ✓ VALID - Object meets all requirements")
        else:
            print("RESULT: ✗ INVALID - Object has validation errors")
        print("="*60 + "\n")


def find_seurat_in_list(obj, parent_name: str = "object") -> Tuple[any, str, str]:
    """
    Recursively search for Seurat object in a list

    Returns:
        Tuple of (Seurat object, path string, type string)
    """
    obj_class = str(ro.r['class'](obj)[0])

    # If this is a Seurat object, return it
    if obj_class == "Seurat":
        return obj, parent_name, obj_class

    # If this is a list, search inside
    if obj_class == "list":
        # Try to get names
        try:
            names = list(ro.r['names'](obj))
        except:
            names = None

        # Search each element
        for i in range(len(obj)):
            if names and i < len(names):
                element_name = names[i]
                full_path = f"{parent_name}${element_name}"
            else:
                element_name = None
                full_path = f"{parent_name}[[{i+1}]]"

            try:
                element = obj[i]
                element_class = str(ro.r['class'](element)[0])

                if element_class == "Seurat":
                    return element, full_path, element_class
                elif element_class == "list":
                    # Recursively search
                    result = find_seurat_in_list(element, full_path)
                    if result[0] is not None:
                        return result
            except:
                continue

    # Not found
    return None, None, None


def load_r_object(file_path: Path, object_name: Optional[str] = None) -> Tuple[any, str]:
    """
    Load R object from .rds or .Rdata file
    Automatically searches for Seurat objects inside lists.

    Returns:
        Tuple of (R object, object type string)
    """
    file_path = Path(file_path).expanduser().resolve()

    if not file_path.exists():
        raise FileNotFoundError(f"File not found: {file_path}")

    suffix = file_path.suffix.lower()

    if suffix == '.rds':
        obj = ro.r['readRDS'](str(file_path))
        obj_type = str(ro.r['class'](obj)[0])

        # If it's a list, try to find Seurat object inside
        if obj_type == "list":
            seurat_obj, path, seurat_type = find_seurat_in_list(obj, "RDS")
            if seurat_obj is not None:
                print(f"Found Seurat object at: {path}")
                return seurat_obj, seurat_type
            else:
                raise ValueError(
                    f"RDS file contains a list but no Seurat object found inside"
                )

        return obj, obj_type

    elif suffix in ['.rdata', '.rda']:
        env = ro.r['new.env']()
        ro.r['load'](str(file_path), env)
        obj_names = list(ro.r['ls'](env))

        if len(obj_names) == 0:
            raise ValueError(f"No objects found in {file_path}")

        # Select object
        if len(obj_names) > 1:
            if object_name is None:
                raise ValueError(
                    f"Multiple objects found in {file_path}: {obj_names}. "
                    f"Please specify --object-name"
                )
            if object_name not in obj_names:
                raise ValueError(
                    f"Object '{object_name}' not found. Available: {obj_names}"
                )
            selected_obj_name = object_name
        else:
            selected_obj_name = obj_names[0]
            if object_name and object_name != obj_names[0]:
                raise ValueError(
                    f"Object '{object_name}' not found. Found: {obj_names[0]}"
                )

        obj = env[selected_obj_name]
        obj_type = str(ro.r['class'](obj)[0])

        # If it's a list, try to find Seurat object inside
        if obj_type == "list":
            seurat_obj, path, seurat_type = find_seurat_in_list(obj, selected_obj_name)
            if seurat_obj is not None:
                print(f"Found Seurat object at: {path}")
                return seurat_obj, seurat_type
            else:
                raise ValueError(
                    f"Object '{selected_obj_name}' is a list but no Seurat object found inside"
                )

        return obj, obj_type

    else:
        raise ValueError(f"Unsupported file format: {suffix}. Use .rds or .Rdata")


def detect_data_type(obj) -> str:
    """
    Automatically detect whether the object is scRNA or scATAC based on assay structure.

    Detection logic:
    - If has 'ACTIVITY' assay → scATAC
    - If has 'RNA', 'SCT', or similar assays with integer counts → scRNA
    - Otherwise, default to scRNA

    Returns:
        'scrna' or 'scatac'
    """
    try:
        # Get all assay names
        assay_names = list(ro.r['names'](ro.r('function(obj) obj@assays')(obj)))

        # Check for ACTIVITY assay (scATAC marker)
        if 'ACTIVITY' in assay_names:
            return 'scatac'

        # Check for typical scRNA assays
        scrna_assays = ['RNA', 'SCT', 'integrated']
        if any(assay in assay_names for assay in scrna_assays):
            return 'scrna'

        # Default to scRNA if uncertain
        return 'scrna'

    except Exception as e:
        # Default to scRNA on any error
        print(f"Warning: Could not auto-detect data type, assuming scRNA ({e})")
        return 'scrna'


def sample_matrix(matrix, max_samples: int = 10000) -> np.ndarray:
    """Sample matrix for efficient checking"""
    try:
        mat_shape = tuple(ro.r['dim'](matrix))
        total_elements = mat_shape[0] * mat_shape[1]

        if total_elements <= max_samples:
            # Small enough, convert all
            return np.array(ro.r['as.matrix'](matrix))
        else:
            # Sample random cells and genes using R indexing
            n_sample_cells = min(100, mat_shape[1])
            n_sample_genes = min(100, mat_shape[0])

            cell_indices = np.random.choice(mat_shape[1], n_sample_cells, replace=False) + 1  # R is 1-indexed
            gene_indices = np.random.choice(mat_shape[0], n_sample_genes, replace=False) + 1

            # Use R function to subset the matrix
            subset_func = ro.r('function(m, rows, cols) as.matrix(m[rows, cols])')
            subset = subset_func(matrix,
                               ro.IntVector(gene_indices.tolist()),
                               ro.IntVector(cell_indices.tolist()))
            return np.array(subset)
    except Exception as e:
        # Fallback: try direct conversion
        try:
            return np.array(ro.r['as.matrix'](matrix))
        except:
            raise Exception(f"Cannot sample matrix: {str(e)}")


def check_matrix_properties(matrix, sample_size: int = 10000) -> Dict[str, any]:
    """Check matrix properties: non-negative, integer ratio, value range"""
    try:
        sampled = sample_matrix(matrix, sample_size).flatten()
        sampled = sampled[~np.isnan(sampled)]  # Remove NaNs

        if len(sampled) == 0:
            return {"error": "Matrix is empty or all NaN"}

        min_val = float(np.min(sampled))
        max_val = float(np.max(sampled))
        mean_val = float(np.mean(sampled))

        # Check if values are non-negative
        non_negative = min_val >= 0

        # Calculate integer ratio
        int_ratio = float(np.mean(np.abs(sampled - np.round(sampled)) < 1e-6))

        return {
            "min": min_val,
            "max": max_val,
            "mean": mean_val,
            "int_ratio": int_ratio,
            "non_negative": non_negative,
            "n_samples": len(sampled)
        }
    except Exception as e:
        return {"error": str(e)}


def get_counts_from_assay(assay):
    """
    Get counts matrix from assay, supporting both Seurat v3/v4 and v5 structures

    Returns matrix or None if not found
    """
    # Check if this is Seurat v5 (has layers slot)
    try:
        slot_names = list(ro.r['slotNames'](assay))

        if 'layers' in slot_names:
            # Seurat v5: get counts from layers
            layers = ro.r('function(a) slot(a, "layers")')(assay)
            layer_names = list(ro.r['names'](layers))

            if 'counts' in layer_names:
                # Access the counts layer
                counts = ro.r('function(l) l[["counts"]]')(layers)
                if int(ro.r['length'](counts)[0]) > 0:
                    return counts

        # Seurat v3/v4: get counts from counts slot
        if 'counts' in slot_names:
            counts = ro.r('function(a) slot(a, "counts")')(assay)
            if int(ro.r['length'](counts)[0]) > 0:
                return counts
    except:
        pass

    return None


def validate_scrna_object(obj, result: ValidationResult, sample_size: int = 10000):
    """Validate scRNA Seurat object"""

    # Check object type
    obj_class = str(ro.r['class'](obj)[0])
    result.info["Object Type"] = obj_class

    if obj_class != "Seurat":
        result.add_error(f"Object is not a Seurat object (found: {obj_class})")
        return

    # Check Seurat version
    seurat_version = "unknown"
    try:
        version_info = ro.r('function(obj) slot(obj, "version")')(obj)
        # Convert to character to get proper version string (e.g., "5.0.2")
        version = str(ro.r['as.character'](version_info)[0])
        result.info["Seurat Version"] = version
        seurat_version = version

        if version.startswith(('3.', '4.', '5.')):
            result.add_pass(f"Seurat version {version} is supported (v3-v5)")
        else:
            result.add_warning(f"Seurat version {version} may not be supported (expected v3-v5)")
    except:
        result.add_warning("Cannot determine Seurat version")

    # Get assay names
    assay_names = list(ro.r['names'](ro.r['slot'](obj, 'assays')))
    result.info["Assays"] = ", ".join(assay_names)

    if not assay_names:
        result.add_error("No assays found in object")
        return

    # Check for counts matrix
    counts_assays = []
    rejected_assays = []
    failed_assays = []
    is_seurat_v5 = seurat_version.startswith('5.')

    for assay_name in assay_names:
        try:
            # Get assay using R function
            assay = ro.r(f'function(obj) obj@assays[["{assay_name}"]]')(obj)

            # Get counts matrix
            counts = get_counts_from_assay(assay)

            if counts is not None:
                # Check matrix properties
                props = check_matrix_properties(counts, sample_size)

                if "error" not in props:
                    if props["non_negative"] and props["int_ratio"] >= 0.98:
                        counts_assays.append({
                            "name": assay_name,
                            "int_ratio": props["int_ratio"],
                            "min": props["min"],
                            "max": props["max"],
                            "mean": props["mean"]
                        })
                    else:
                        # Record why it was rejected
                        reasons = []
                        if not props["non_negative"]:
                            reasons.append(f"has negative values (min={props['min']:.4f})")
                        if props["int_ratio"] < 0.98:
                            reasons.append(f"int_ratio={props['int_ratio']:.3f} < 0.98 (appears normalized)")
                        rejected_assays.append({
                            "name": assay_name,
                            "reasons": reasons
                        })
                else:
                    failed_assays.append({
                        "name": assay_name,
                        "error": props["error"]
                    })
            else:
                failed_assays.append({
                    "name": assay_name,
                    "error": "counts matrix not found"
                })
        except Exception as e:
            failed_assays.append({
                "name": assay_name,
                "error": str(e)
            })

    if counts_assays:
        result.add_pass(f"Found {len(counts_assays)} assay(s) with valid counts matrix")
        for assay_info in counts_assays:
            result.add_pass(
                f"  - {assay_info['name']}: int_ratio={assay_info['int_ratio']:.3f}, "
                f"range=[{assay_info['min']:.1f}, {assay_info['max']:.1f}]"
            )
    else:
        error_msg = "No usable counts assay found (need non-negative matrix with int_ratio ≥ 0.98)"
        if rejected_assays:
            error_msg += f"\nRejected {len(rejected_assays)} assay(s):"
            for rej in rejected_assays:
                error_msg += f"\n  - {rej['name']}: {', '.join(rej['reasons'])}"
        if failed_assays:
            error_msg += f"\nFailed to check {len(failed_assays)} assay(s):"
            for fail in failed_assays:
                error_msg += f"\n  - {fail['name']}: {fail['error']}"
        result.add_error(error_msg)

    # Check reductions (UMAP/tSNE)
    try:
        reductions = ro.r['names'](ro.r['slot'](obj, 'reductions'))
        if len(reductions) > 0:
            reduction_list = list(reductions)
            result.info["Reductions"] = ", ".join(reduction_list)

            has_umap = 'umap' in [r.lower() for r in reduction_list]
            has_tsne = 'tsne' in [r.lower() for r in reduction_list]

            if has_umap or has_tsne:
                result.add_pass(f"Dimensionality reduction found: {', '.join(reduction_list)}")
            else:
                result.add_warning(
                    f"No standard UMAP/tSNE reduction found (found: {', '.join(reduction_list)})"
                )

            # Check if embeddings row names match cell names
            for red_name in reduction_list:
                try:
                    reduction = ro.r(f'function(obj) obj@reductions[["{red_name}"]]')(obj)
                    embeddings = ro.r('function(r) slot(r, "cell.embeddings")')(reduction)

                    # Get embedding rownames
                    emb_rownames_result = ro.r['rownames'](embeddings)
                    if ro.r['is.null'](emb_rownames_result)[0]:
                        result.add_warning(f"  - {red_name}: embeddings have no row names")
                        continue
                    emb_rownames = list(emb_rownames_result)

                    # Get cell names from Seurat object
                    # For Seurat v5, use Cells() function; for v3/v4, get from matrix
                    try:
                        # Use Cells() function (works for all versions)
                        cell_names = list(ro.r('Cells')(obj))
                    except:
                        # Fallback: get from matrix colnames
                        first_assay = ro.r(f'function(obj) obj@assays[["{assay_names[0]}"]]')(obj)
                        counts = get_counts_from_assay(first_assay)
                        if counts is not None:
                            colnames_result = ro.r['colnames'](counts)
                            if not ro.r['is.null'](colnames_result)[0]:
                                cell_names = list(colnames_result)
                            else:
                                result.add_warning(f"  - {red_name}: cannot get cell names")
                                continue
                        else:
                            result.add_warning(f"  - {red_name}: cannot get cell names")
                            continue

                    if emb_rownames == cell_names:
                        result.add_pass(f"  - {red_name}: embeddings row names match cell names")
                    else:
                        result.add_error(
                            f"  - {red_name}: embeddings row names DO NOT match cell names "
                            f"(found {len(emb_rownames)}, expected {len(cell_names)})"
                        )
                except Exception as e:
                    result.add_warning(f"  - {red_name}: cannot verify row names ({str(e)})")
        else:
            result.add_warning("No dimensionality reductions found (UMAP/tSNE recommended)")
    except Exception as e:
        result.add_warning(f"Cannot check reductions: {str(e)}")

    # Check meta.data
    try:
        meta_data = ro.r('function(obj) slot(obj, "meta.data")')(obj)
        n_cells_meta = ro.r['nrow'](meta_data)[0]
        meta_cols = list(ro.r['colnames'](meta_data))

        result.info["N cells (meta.data)"] = n_cells_meta
        result.info["Meta columns"] = len(meta_cols)

        # Check for clustering column
        cluster_cols = [col for col in meta_cols if 'cluster' in col.lower() or 'seurat_clusters' in col.lower()]
        if cluster_cols:
            result.add_pass(f"Clustering column(s) found: {', '.join(cluster_cols)}")
        else:
            result.add_warning("No obvious clustering column found (e.g., 'seurat_clusters')")

        # Check meta.data row names match cell names
        meta_rownames_result = ro.r['rownames'](meta_data)
        if not ro.r['is.null'](meta_rownames_result)[0]:
            meta_rownames = list(meta_rownames_result)

            # Get cell names using Cells() function
            try:
                cell_names = list(ro.r('Cells')(obj))
            except:
                # Fallback: try to get from counts matrix
                first_assay = ro.r(f'function(obj) obj@assays[["{assay_names[0]}"]]')(obj)
                counts = get_counts_from_assay(first_assay)
                if counts is not None:
                    colnames_result = ro.r['colnames'](counts)
                    if not ro.r['is.null'](colnames_result)[0]:
                        cell_names = list(colnames_result)
                    else:
                        result.add_warning("Cannot get cell names for meta.data comparison")
                        cell_names = None
                else:
                    result.add_warning("Cannot get cell names for meta.data comparison")
                    cell_names = None

            if cell_names is not None:
                if meta_rownames == cell_names:
                    result.add_pass("meta.data row names match cell names")
                else:
                    result.add_error(
                        f"meta.data row names DO NOT match cell names "
                        f"(found {len(meta_rownames)}, expected {len(cell_names)})"
                    )
        else:
            result.add_warning("meta.data has no row names")
    except Exception as e:
        result.add_error(f"Cannot check meta.data: {str(e)}")

    # Check gene annotations (meta.features)
    try:
        first_assay = ro.r(f'function(obj) obj@assays[["{assay_names[0]}"]]')(obj)

        # For Seurat v5, meta.data is used instead of meta.features
        if is_seurat_v5:
            try:
                meta_features = ro.r('function(a) slot(a, "meta.data")')(first_assay)
                n_features = ro.r['nrow'](meta_features)[0]
                feature_cols = list(ro.r['colnames'](meta_features))

                if n_features > 0 and len(feature_cols) > 0:
                    result.add_pass(f"Gene metadata (meta.data) found with {len(feature_cols)} columns")
                    if 'gene_name' in feature_cols:
                        result.add_pass("  - 'gene_name' column present")
                else:
                    result.add_warning("Gene metadata is empty (gene annotations recommended)")
            except:
                result.add_warning("Cannot access gene metadata (annotations recommended)")
        else:
            # Seurat v3/v4
            try:
                meta_features = ro.r('function(a) slot(a, "meta.features")')(first_assay)
                n_features = ro.r['nrow'](meta_features)[0]
                feature_cols = list(ro.r['colnames'](meta_features))

                if n_features > 0 and len(feature_cols) > 0:
                    result.add_pass(f"meta.features found with {len(feature_cols)} columns")
                    if 'gene_name' in feature_cols:
                        result.add_pass("  - 'gene_name' column present")
                else:
                    result.add_warning("meta.features is empty (gene annotations recommended)")
            except:
                result.add_warning("Cannot access meta.features (gene annotations recommended)")
    except Exception as e:
        result.add_warning(f"Cannot check gene annotations: {str(e)}")

    # Check gene ID format (row names)
    try:
        # Get gene names using Features() function (works for all versions)
        try:
            gene_names = list(ro.r('Features')(obj, assay=assay_names[0]))
        except:
            # Fallback: try to get from counts matrix
            first_assay = ro.r(f'function(obj) obj@assays[["{assay_names[0]}"]]')(obj)
            counts = get_counts_from_assay(first_assay)
            if counts is not None:
                rownames_result = ro.r['rownames'](counts)
                if not ro.r['is.null'](rownames_result)[0]:
                    gene_names = list(rownames_result)
                else:
                    result.add_warning("Cannot get gene names from counts matrix")
                    gene_names = None
            else:
                result.add_warning("Cannot get gene names from counts matrix")
                gene_names = None

        if gene_names is not None:
            result.info["N genes"] = len(gene_names)

            # Check for Ensembl IDs
            ensembl_pattern = ['ENSDARG', 'ENSMUSG', 'ENSG']  # zebrafish, mouse, human
            has_ensembl = any(any(gene.startswith(pattern) for pattern in ensembl_pattern)
                             for gene in gene_names[:100])  # Check first 100

            if has_ensembl:
                result.add_pass("Gene names appear to use Ensembl IDs")
            else:
                result.add_warning(
                    "Gene names may not use standard Ensembl IDs "
                    "(ensure consistency with species adapter)"
                )
    except Exception as e:
        result.add_warning(f"Cannot check gene names: {str(e)}")


def get_data_from_assay(obj, assay_name: str, slot_name: str = "data"):
    """
    Get data matrix from assay, supporting both Seurat v3/v4 and v5 structures

    Returns matrix or None if not found
    """
    try:
        # Get assay
        assay = ro.r(f'function(obj) obj@assays[["{assay_name}"]]')(obj)

        # Check if this is Seurat v5 (has layers slot)
        slot_names = list(ro.r['slotNames'](assay))

        if 'layers' in slot_names:
            # Seurat v5: get data from layers
            layers = ro.r('function(a) slot(a, "layers")')(assay)
            layer_names = list(ro.r['names'](layers))

            if slot_name in layer_names:
                # Access the layer
                data = ro.r(f'function(l) l[["{slot_name}"]]')(layers)
                if int(ro.r['length'](data)[0]) > 0:
                    return data

        # Seurat v3/v4: get from slot directly
        if slot_name in slot_names:
            data = ro.r(f'function(a) slot(a, "{slot_name}")')(assay)
            if int(ro.r['length'](data)[0]) > 0:
                return data
    except:
        pass

    return None


def validate_scatac_object(obj, result: ValidationResult,
                          assay_name: str = "ACTIVITY",
                          slot_name: str = "data",
                          sample_size: int = 10000):
    """Validate scATAC Seurat/Signac object"""

    # Check object type
    obj_class = str(ro.r['class'](obj)[0])
    result.info["Object Type"] = obj_class

    if obj_class != "Seurat":
        result.add_error(f"Object is not a Seurat object (found: {obj_class})")
        return

    # Get assay names
    assay_names = list(ro.r['names'](ro.r['slot'](obj, 'assays')))
    result.info["Assays"] = ", ".join(assay_names)

    if not assay_names:
        result.add_error("No assays found in object")
        return

    # Check for specified assay
    if assay_name not in assay_names:
        result.add_error(f"Assay '{assay_name}' not found (available: {', '.join(assay_names)})")
        return
    else:
        result.add_pass(f"Assay '{assay_name}' found")

    # Check for specified slot/layer
    try:
        data_matrix = get_data_from_assay(obj, assay_name, slot_name)

        if data_matrix is None:
            result.add_error(f"Layer/slot '{slot_name}' in assay '{assay_name}' not found or empty")
            return
        else:
            result.add_pass(f"Layer/slot '{slot_name}' found in assay '{assay_name}'")

        # Check matrix properties
        props = check_matrix_properties(data_matrix, sample_size)

        if "error" in props:
            result.add_error(f"Cannot check matrix properties: {props['error']}")
        else:
            result.info["Data min"] = f"{props['min']:.4f}"
            result.info["Data max"] = f"{props['max']:.4f}"
            result.info["Data mean"] = f"{props['mean']:.4f}"
            result.info["Integer ratio"] = f"{props['int_ratio']:.4f}"

            # Validate: should be normalized (not counts)
            if props["int_ratio"] > 0.9:
                result.add_error(
                    f"Data appears to be counts (int_ratio={props['int_ratio']:.3f}). "
                    "Expected log1p normalized Gene Activity."
                )
            else:
                result.add_pass(f"Data appears normalized (int_ratio={props['int_ratio']:.3f})")

            if props["max"] > 50:
                result.add_error(
                    f"Data max value is {props['max']:.2f} (should be < 50 for normalized data)"
                )
            else:
                result.add_pass(f"Data max value is within range ({props['max']:.2f} < 50)")

            if props["min"] < 0:
                result.add_error(
                    f"Data contains negative values (min={props['min']:.4f}). "
                    "This suggests scaled data, not log1p normalized."
                )
            else:
                result.add_pass("Data is non-negative")

    except Exception as e:
        result.add_error(f"Cannot access layer/slot '{slot_name}': {str(e)}")
        return

    # Check for UMAP reduction (required for scATAC)
    try:
        reductions = ro.r['names'](ro.r['slot'](obj, 'reductions'))
        if len(reductions) > 0:
            reduction_list = list(reductions)
            result.info["Reductions"] = ", ".join(reduction_list)

            # Look for any UMAP reduction (case-insensitive, partial match)
            umap_reductions = [r for r in reduction_list if 'umap' in r.lower()]

            if umap_reductions:
                umap_name = umap_reductions[0]  # Use the first UMAP found
                result.add_pass(f"UMAP reduction found: '{umap_name}'")

                # Check if embeddings row names match cell names
                try:
                    reduction = ro.r(f'function(obj) obj@reductions[["{umap_name}"]]')(obj)
                    embeddings = ro.r('function(r) slot(r, "cell.embeddings")')(reduction)

                    # Get embedding rownames
                    emb_rownames_result = ro.r['rownames'](embeddings)
                    if ro.r['is.null'](emb_rownames_result)[0]:
                        result.add_warning(f"UMAP embeddings have no row names")
                    else:
                        emb_rownames = list(emb_rownames_result)

                        # Get cell names using Cells() function
                        try:
                            cell_names = list(ro.r('Cells')(obj))
                        except:
                            # Fallback: try to get from data matrix
                            colnames_result = ro.r['colnames'](data_matrix)
                            if not ro.r['is.null'](colnames_result)[0]:
                                cell_names = list(colnames_result)
                            else:
                                result.add_warning("Cannot get cell names for UMAP comparison")
                                cell_names = None

                        if cell_names is not None:
                            if emb_rownames == cell_names:
                                result.add_pass("UMAP embeddings row names match cell names")
                            else:
                                result.add_error(
                                    f"UMAP embeddings row names DO NOT match cell names "
                                    f"(found {len(emb_rownames)}, expected {len(cell_names)})"
                                )
                except Exception as e:
                    result.add_warning(f"Cannot verify UMAP row names: {str(e)}")
            else:
                result.add_error("UMAP reduction not found (required for scATAC)")
        else:
            result.add_error("No reductions found (UMAP required for scATAC)")
    except Exception as e:
        result.add_error(f"Cannot check reductions: {str(e)}")

    # Check meta.data
    try:
        meta_data = ro.r('function(obj) slot(obj, "meta.data")')(obj)
        n_cells_meta = ro.r['nrow'](meta_data)[0]
        meta_cols = list(ro.r['colnames'](meta_data))

        result.info["N cells (meta.data)"] = n_cells_meta
        result.info["Meta columns"] = len(meta_cols)

        # Check for grouping columns
        grouping_cols = [col for col in meta_cols
                        if any(keyword in col.lower()
                              for keyword in ['cluster', 'cell_type', 'stage', 'condition'])]
        if grouping_cols:
            result.add_pass(f"Grouping column(s) found: {', '.join(grouping_cols)}")
        else:
            result.add_warning("No obvious grouping columns found (recommended for visualization)")

        # Check meta.data row names match cell names
        meta_rownames_result = ro.r['rownames'](meta_data)
        if not ro.r['is.null'](meta_rownames_result)[0]:
            meta_rownames = list(meta_rownames_result)

            # Get cell names using Cells() function
            try:
                cell_names = list(ro.r('Cells')(obj))
            except:
                # Fallback: try to get from data matrix
                colnames_result = ro.r['colnames'](data_matrix)
                if not ro.r['is.null'](colnames_result)[0]:
                    cell_names = list(colnames_result)
                else:
                    result.add_warning("Cannot get cell names for meta.data comparison")
                    cell_names = None

            if cell_names is not None:
                if meta_rownames == cell_names:
                    result.add_pass("meta.data row names match cell names")
                else:
                    result.add_error(
                        f"meta.data row names DO NOT match cell names "
                        f"(found {len(meta_rownames)}, expected {len(cell_names)})"
                    )
        else:
            result.add_warning("meta.data has no row names")
    except Exception as e:
        result.add_error(f"Cannot check meta.data: {str(e)}")

    # Check gene annotations (meta.features)
    try:
        # Get gene names using Features() function
        try:
            gene_names = list(ro.r('Features')(obj, assay=assay_name))
        except:
            # Fallback: try to get from data matrix rownames
            rownames_result = ro.r['rownames'](data_matrix)
            if not ro.r['is.null'](rownames_result)[0]:
                gene_names = list(rownames_result)
            else:
                gene_names = None

        if gene_names is not None:
            result.info["N genes"] = len(gene_names)
            if len(gene_names) > 0:
                result.add_pass("Gene IDs found")
            else:
                result.add_warning("No gene IDs found")
        else:
            result.add_warning("Cannot get gene IDs")
    except Exception as e:
        result.add_warning(f"Cannot check gene annotations: {str(e)}")


def main():
    parser = argparse.ArgumentParser(
        description="Validate Seurat objects against scRNA/scATAC conversion requirements",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Auto-detect data type
  python validate_seurat_object.py --input sample.rds

  # Validate scRNA object (explicit)
  python validate_seurat_object.py --input sample.rds --type scrna

  # Validate scATAC object (explicit)
  python validate_seurat_object.py --input sample.rds --type scatac

  # Validate scATAC with custom assay/slot
  python validate_seurat_object.py --input sample.rds --type scatac \\
      --assay ChromatinAssay --slot data

  # Validate RData with multiple objects
  python validate_seurat_object.py --input sample.RData --type scrna \\
      --object-name seurat_obj
        """
    )

    parser.add_argument(
        '--input', '-i',
        required=True,
        help="Path to .rds or .RData file"
    )

    parser.add_argument(
        '--type', '-t',
        choices=['scrna', 'scatac'],
        default=None,
        help="Data type: scrna or scatac (auto-detect if not specified)"
    )

    parser.add_argument(
        '--object-name',
        help="Object name (required for .RData with multiple objects)"
    )

    parser.add_argument(
        '--assay',
        default='ACTIVITY',
        help="Assay name for scATAC (default: ACTIVITY)"
    )

    parser.add_argument(
        '--slot',
        default='data',
        help="Slot name for scATAC (default: data)"
    )

    parser.add_argument(
        '--sample-size',
        type=int,
        default=10000,
        help="Number of matrix elements to sample for checking (default: 10000)"
    )

    parser.add_argument(
        '--verbose', '-v',
        action='store_true',
        help="Show all passed checks"
    )

    parser.add_argument(
        '--show-r-output',
        action='store_true',
        help="Show R package loading messages (hidden by default)"
    )

    parser.add_argument(
        '--json',
        help="Output validation results to JSON file"
    )

    args = parser.parse_args()

    # Setup R environment (suppress output unless --show-r-output is specified)
    setup_r_environment(quiet=not args.show_r_output)

    # Initialize result
    result = ValidationResult()

    try:
        # Load Seurat library in R to enable Cells() and Features() functions
        try:
            # Use suppressPackageStartupMessages to hide loading messages
            ro.r('suppressPackageStartupMessages(library(Seurat))')
        except:
            pass  # Seurat may already be loaded or not available

        # Load object
        print(f"Loading object from {args.input}...")
        obj, obj_type = load_r_object(args.input, args.object_name)

        # Auto-detect data type if not specified
        if args.type is None:
            data_type = detect_data_type(obj)
            print(f"Auto-detected data type: {data_type.upper()}")
        else:
            data_type = args.type

        result.info["File"] = str(args.input)
        result.info["Data Type"] = data_type.upper()

        # Validate based on type
        if data_type == 'scrna':
            validate_scrna_object(obj, result, args.sample_size)
        elif data_type == 'scatac':
            validate_scatac_object(obj, result, args.assay, args.slot, args.sample_size)

    except Exception as e:
        result.add_error(f"Fatal error: {str(e)}")
        import traceback
        if args.verbose:
            traceback.print_exc()

    # Print summary
    result.print_summary(verbose=args.verbose)

    # Save JSON if requested
    if args.json:
        output_data = {
            "valid": result.is_valid(),
            "info": result.info,
            "passed": result.passed,
            "warnings": result.warnings,
            "errors": result.errors
        }
        with open(args.json, 'w') as f:
            json.dump(output_data, f, indent=2)
        print(f"Results saved to {args.json}")

    # Exit code
    sys.exit(0 if result.is_valid() else 1)


if __name__ == '__main__':
    main()
