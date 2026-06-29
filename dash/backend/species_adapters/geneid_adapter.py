"""
Adapter for species that use a generic GeneID-style primary identifier.

This keeps the StandardAdapter query logic while fixing user-facing labels
that would otherwise incorrectly say "Ensembl ID".
"""

from __future__ import annotations

from typing import Any, Dict, List

from .standard_adapter import StandardAdapter


class GeneIdAdapter(StandardAdapter):
    """Adapter for GeneID-based species such as killifish."""

    def get_available_search_fields(self) -> List[Dict[str, str]]:
        fields = super().get_available_search_fields()
        return [
            {"label": "Gene ID", "value": field["value"]}
            if field.get("value") == "ens_id"
            else field
            for field in fields
        ]

    def get_search_fields(self) -> List[Dict[str, str]]:
        return [
            {"label": "Gene Name", "value": "gene_name"},
            {"label": "Gene ID", "value": "ens_id"},
            {"label": "Description", "value": "description"},
            {"label": "GO Accession", "value": "go_id"},
            {"label": "GO Term", "value": "go_term"},
        ]

    def format_gene_label(self, row) -> str:
        primary_id = row.get("ENS", "")
        gene_name = row.get("GeneName", "")
        if not gene_name or gene_name == primary_id:
            return primary_id
        return f"{gene_name} ({primary_id})"

    def get_external_links(self, result: Dict[str, Any]) -> Dict[str, str]:
        primary_id = result.get("primary_id")
        if primary_id and primary_id.isdigit():
            return {"NCBI Gene": f"https://www.ncbi.nlm.nih.gov/gene/{primary_id}"}
        return {}

    def get_tf_families(self, db_conn, primary_id: str, secondary_id: str = None) -> List[str]:
        cursor = db_conn.cursor()
        table_exists = cursor.execute(
            "SELECT 1 FROM sqlite_master WHERE type = 'table' AND name = 'ENS_TF'"
        ).fetchone()
        if not table_exists:
            return []
        return super().get_tf_families(db_conn, primary_id, secondary_id)

    def get_summary_columns(self) -> Dict[str, str]:
        return {
            "gene_name": "Gene Name",
            "primary_id": "Gene ID",
            "locus": "Locus",
            "tf_families": "TF Family",
            "description": "Description",
        }

    def get_column_display_config(self) -> Dict[str, str]:
        return {
            "gene_name": "Gene Name",
            "primary_id": "Gene ID",
            "secondary_id": "Secondary ID",
            "locus": "Locus",
            "description": "Description",
            "gene_type": "Gene Type",
            "tf_families": "TF Family",
        }

    def format_search_result_display(self, result: Dict[str, Any]) -> Dict[str, Any]:
        gene_name = result.get("gene_name", "")
        primary_id = result.get("primary_id", "")

        if not gene_name or gene_name == primary_id:
            return {
                "main_text": primary_id or gene_name or "N/A",
                "secondary_items": [],
            }

        return {
            "main_text": gene_name,
            "secondary_items": [
                {"text": primary_id, "className": "text-muted"} if primary_id else None
            ],
        }
