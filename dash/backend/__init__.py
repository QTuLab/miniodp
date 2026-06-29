"""
Backend API module for Multi-Omics Data Portal
Provides data access interfaces for the Dash application
"""

from .gene_search import GeneSearchAPI

__all__ = ['GeneSearchAPI']