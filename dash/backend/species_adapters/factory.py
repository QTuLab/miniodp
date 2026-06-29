"""
Configuration-driven factory for creating species adapters.

This module provides a zero species-specific code factory that uses
configuration files to determine adapter selection and behavior.
"""
import logging
from pathlib import Path
from typing import Dict, Any, Optional

try:
    import tomllib
except ModuleNotFoundError:  # Python < 3.11
    import tomli as tomllib

from .base_adapter import SpeciesAdapter
from .standard_adapter import StandardAdapter
from .medaka_adapter import MedakaAdapter
from .geneid_adapter import GeneIdAdapter


logger = logging.getLogger(__name__)


# Adapter registry - maps adapter type names to classes
ADAPTER_REGISTRY = {
    'standard': StandardAdapter,
    'dual_system': MedakaAdapter,
    'gene_id': GeneIdAdapter,
}


def _load_species_configuration() -> Dict[str, Any]:
    """
    Load species configuration from TOML file.
    
    Returns:
        Configuration dictionary
    """
    config_path = Path(__file__).parent.parent.parent / "data" / "species_adapters.toml"
    
    try:
        with open(config_path, "rb") as f:
            return tomllib.load(f)
    except Exception as e:
        logger.warning("Could not load species configuration: %s", e)
        # Return minimal fallback configuration
        return {
            'default': {'adapter_type': 'standard'},
            'adapters': {'standard': 'StandardAdapter'},
            'species': {}
        }


def create_species_adapter(species_key: str) -> SpeciesAdapter:
    """
    Create the appropriate species adapter based on configuration.
    
    This is a zero species-specific code implementation that uses
    configuration to determine adapter selection.
    
    Args:
        species_key: Identifier for the species
        
    Returns:
        Appropriate SpeciesAdapter instance
        
    Raises:
        ValueError: If adapter type is not found in registry
    """
    if not species_key:
        raise ValueError("Species key cannot be empty")
    
    species_key = species_key.lower().strip()
    config = _load_species_configuration()
    
    # Start with default configuration and merge with species-specific overrides
    default_config = config.get('default', {})
    species_specific = config.get('species', {}).get(species_key, {})
    
    # Merge configurations (species-specific overrides defaults)
    merged_config = default_config.copy()
    merged_config.update(species_specific)
    
    adapter_type = merged_config.get('adapter_type', 'standard')
    
    # Look up adapter class in registry
    if adapter_type not in ADAPTER_REGISTRY:
        available_types = list(ADAPTER_REGISTRY.keys())
        raise ValueError(f"Unknown adapter type '{adapter_type}'. Available: {available_types}")
    
    adapter_class = ADAPTER_REGISTRY[adapter_type]
    
    # Create adapter instance with configuration
    return adapter_class(species_key)


def get_supported_species() -> Dict[str, str]:
    """
    Get list of all configured species and their adapter types.
    
    Returns:
        Dict mapping species keys to adapter type descriptions
    """
    config = _load_species_configuration()
    species_config = config.get('species', {})
    
    result = {}
    for species_key, species_data in species_config.items():
        adapter_type = species_data.get('adapter_type', 'standard')
        adapter_class_name = config.get('adapters', {}).get(adapter_type, 'Unknown')
        result[species_key] = adapter_class_name
    
    return result


def get_species_configuration(species_key: str) -> Optional[Dict[str, Any]]:
    """
    Get configuration for a specific species.
    
    Args:
        species_key: Species identifier
        
    Returns:
        Species configuration dict or None if not found
    """
    config = _load_species_configuration()
    return config.get('species', {}).get(species_key.lower().strip())


def validate_species_key(species_key: str) -> bool:
    """
    Check if a species key is configured.
    
    Args:
        species_key: Species identifier to validate
        
    Returns:
        True if configured, False otherwise
    """
    config = _load_species_configuration()
    return species_key.lower().strip() in config.get('species', {})


def add_species_configuration(species_key: str, adapter_type: str, **kwargs) -> bool:
    """
    Programmatically add a new species configuration.
    
    This allows runtime addition of species without modifying the TOML file.
    Note: This only affects the current session.
    
    Args:
        species_key: Species identifier
        adapter_type: Adapter type to use
        **kwargs: Additional configuration parameters
        
    Returns:
        True if successfully added, False otherwise
    """
    if adapter_type not in ADAPTER_REGISTRY:
        return False
    
    # This would require implementing runtime configuration modification
    # For now, return False to indicate this is not yet implemented
    return False


def get_available_adapter_types() -> Dict[str, str]:
    """
    Get list of available adapter types.
    
    Returns:
        Dict mapping adapter type names to class names
    """
    config = _load_species_configuration()
    return config.get('adapters', {})
