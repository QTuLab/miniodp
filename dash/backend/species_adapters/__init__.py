"""
Species-specific adapters for handling different data models and search strategies.

This module provides a clean abstraction layer to handle differences between species
while keeping the core application logic generic and maintainable.
"""

from .factory import create_species_adapter, get_supported_species, validate_species_key
from .base_adapter import SpeciesAdapter
from .standard_adapter import StandardAdapter
from .medaka_adapter import MedakaAdapter
from .geneid_adapter import GeneIdAdapter

__all__ = [
    'create_species_adapter',
    'get_supported_species',
    'validate_species_key',
    'SpeciesAdapter',
    'StandardAdapter',
    'MedakaAdapter',
    'GeneIdAdapter',
]
