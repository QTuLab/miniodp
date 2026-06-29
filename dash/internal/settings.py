import os
from pathlib import Path
from typing import Optional
from urllib.parse import parse_qs, urlencode, urlsplit, urlunsplit

import logging


logger = logging.getLogger("miniodp.dash")


# --- Configuration from Environment Variables ---
# Read configuration from environment variables with fallback to defaults
DASH_PREFIX = os.environ.get("DASH_PREFIX", "/")
DASH_DATA_PATH = os.environ.get("DASH_DATA_PATH", "./data")


# Ensure URL prefix format is correct
if not DASH_PREFIX.startswith("/"):
    DASH_PREFIX = "/" + DASH_PREFIX
if not DASH_PREFIX.endswith("/"):
    DASH_PREFIX = DASH_PREFIX + "/"


def strip_dash_prefix(pathname: Optional[str]) -> Optional[str]:
    """Remove DASH_PREFIX from a pathname if present."""
    if not pathname or DASH_PREFIX == "/":
        return pathname
    normalized_prefix = DASH_PREFIX.rstrip("/")
    if pathname.startswith(normalized_prefix):
        stripped = pathname[len(normalized_prefix) :]
        return stripped or "/"
    return pathname


def build_prefixed_href(relative_path: str = "") -> str:
    """Construct hrefs that always include the configured DASH_PREFIX."""
    stripped = relative_path.lstrip("/")
    return f"{DASH_PREFIX}{stripped}" if stripped else DASH_PREFIX


def build_jbrowse_url(
    base_url: str,
    locus: str,
    padding_mult: int = 2,
    min_span: int = 20000,
    max_span: int = 500000,
) -> Optional[str]:
    """
    Build a JBrowse URL by replacing the loc parameter with a padded locus window.
    """
    if not base_url or not locus:
        return None
    try:
        chrom_part = locus.split(":")[0].replace("chr", "").strip()
        pos_part = locus.split(":")[1].split()[0]  # drop strand info if present
        start_str, end_str = pos_part.replace("-", "..").split("..")
        start = int(start_str)
        end = int(end_str)
        if start > end:
            start, end = end, start

        gene_len = max(1, end - start + 1)
        span = gene_len * (1 + 2 * padding_mult)
        span = max(min_span, min(max_span, span))
        pad = max(0, (span - gene_len) // 2)
        new_start = max(1, start - pad)
        new_end = new_start + span - 1

        split_url = urlsplit(base_url)
        query = parse_qs(split_url.query, keep_blank_values=True)
        query["loc"] = [f"{chrom_part}:{new_start}..{new_end}"]
        new_query = urlencode(query, doseq=True)
        return urlunsplit(
            (split_url.scheme, split_url.netloc, split_url.path, new_query, split_url.fragment)
        )
    except Exception as exc:
        logger.warning(f"Failed to build JBrowse URL from locus {locus}: {exc}")
        return base_url


def load_tf_families_text() -> str:
    """Load TF families markdown/text from assets; returns a human-readable fallback on error."""
    try:
        tf_families_path = Path(__file__).resolve().parent.parent / "assets" / "tf_families.txt"
        return tf_families_path.read_text()
    except Exception as exc:
        logger.error(f"❌ Failed to load TF families: {exc}")
        return "Could not load TF families list."

