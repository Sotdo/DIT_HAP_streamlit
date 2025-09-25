"""
Data configuration module for DIT-HAP Streamlit application.

This module provides centralized configuration for all data file paths
used throughout the application.
"""

# ================================= Imports =================================
from pathlib import Path
from dataclasses import dataclass, field
from typing import Optional, List

# ================================= Constants =================================
POMBASE_DATA_PATH = Path("./data/resource/pombase_data/2025-09-01")
RESOURCE_PATH = Path("./data/resource")
RAW_DATA_PATH = Path("./data/raw/HD_DIT_HAP")
LONG_TIMECOURSE_DATA_PATH = Path("./data/raw/Long_timecourse_data")
HAPLOID_DATA_PATH = Path("./data/raw/haploid_data")

# ================================= Data Configuration Classes =================================

@dataclass
class GeneMetadataConfig:
    """Configuration for gene metadata file paths."""
    
    gene_IDs_names_products: Path = POMBASE_DATA_PATH / "Gene_metadata/gene_IDs_names_products.tsv"
    deletion_library_essentiality: Path = RESOURCE_PATH / "Hayles_2013_OB_merged_categories_sysIDupdated.xlsx"
    genome_intervals: Path = POMBASE_DATA_PATH / "genome_region/genome_intervals.bed"

    def validate_paths(self) -> None:
        """Validate that all required files exist."""
        for field_name, file_path in self.__dict__.items():
            if not file_path.exists():
                raise FileNotFoundError(f"Gene metadata file not found: {file_path}")
            if not file_path.is_file():
                raise ValueError(f"Path is not a file: {file_path}")


@dataclass
class InsertionLevelConfig:
    """Configuration for insertion level statistics file paths."""
    
    insertion_annotations: Path = RAW_DATA_PATH / "insertion_level/annotations.tsv.gz"
    insertion_level_baseMean: Path = RAW_DATA_PATH / "insertion_level/baseMean.tsv"
    imputation_statistics: Optional[Path] = RAW_DATA_PATH / "insertion_level/imputation_statistics.tsv"
    insertion_level_LFCs: Path = RAW_DATA_PATH / "insertion_level/LFC.tsv"
    insertion_level_weights: Path = RAW_DATA_PATH / "insertion_level/transformed_weights.tsv"
    insertion_level_fitting_LFCs: Path = RAW_DATA_PATH / "insertion_level/fitting_LFCs.tsv"
    insertion_level_fitting_results: Path = RAW_DATA_PATH / "insertion_level/fitting_results.tsv"

    def validate_paths(self) -> None:
        """Validate that all required files exist."""
        required_files = [
            self.insertion_annotations,
            self.insertion_level_baseMean,
            self.insertion_level_LFCs,
            self.insertion_level_weights,
            self.insertion_level_fitting_LFCs,
            self.insertion_level_fitting_results
        ]
        
        for file_path in required_files:
            if not file_path.exists():
                raise FileNotFoundError(f"Insertion level file not found: {file_path}")
            if not file_path.is_file():
                raise ValueError(f"Path is not a file: {file_path}")
        
        # Validate optional file if provided
        if self.imputation_statistics and not self.imputation_statistics.exists():
            raise FileNotFoundError(f"Imputation statistics file not found: {self.imputation_statistics}")


@dataclass
class GeneLevelConfig:
    """Configuration for gene level statistics file paths."""
    
    gene_level_LFCs: Path = RAW_DATA_PATH / "gene_level/LFC.tsv"
    gene_level_fitting_LFCs: Path = RAW_DATA_PATH / "gene_level/fitting_LFCs.tsv"
    gene_level_fitting_results: Path = RAW_DATA_PATH / "gene_level/fitting_results.tsv"

    def validate_paths(self) -> None:
        """Validate that all required files exist."""
        for field_name, file_path in self.__dict__.items():
            if not file_path.exists():
                raise FileNotFoundError(f"Gene level file not found: {file_path}")
            if not file_path.is_file():
                raise ValueError(f"Path is not a file: {file_path}")


@dataclass
class GeneOntologyDataConfig:
    """Configuration for ontology data file paths."""
    
    gene_ontology_obo: Path = POMBASE_DATA_PATH / "ontologies_and_associations/go-basic.obo"
    gene_ontology_association_gaf: Path = POMBASE_DATA_PATH / "ontologies_and_associations/gene_ontology_annotation.gaf.tsv"
    goslim_terms_table: List[Path] = field(default_factory=lambda: [
        POMBASE_DATA_PATH / "ontologies_and_associations/bp_go_slim_terms.tsv",
        POMBASE_DATA_PATH / "ontologies_and_associations/cc_go_slim_terms.tsv",
        POMBASE_DATA_PATH / "ontologies_and_associations/mf_go_slim_terms.tsv"
    ])

    def validate_paths(self) -> None:
        """Validate that all required files exist."""
        required_files = [
            self.gene_ontology_obo,
            self.gene_ontology_association_gaf,
            *self.goslim_terms_table
        ]
        
        for file_path in required_files:
            if not file_path.exists():
                raise FileNotFoundError(f"Gene ontology file not found: {file_path}")
            if not file_path.is_file():
                raise ValueError(f"Path is not a file: {file_path}")

@dataclass
class FYPODataConfig:
    """Configuration for fission yeast phenotype ontology data file paths."""
    
    fypo_obo: Path = POMBASE_DATA_PATH / "ontologies_and_associations/fypo-simple-pombase.obo"
    fypo_gaf: Path = POMBASE_DATA_PATH / "ontologies_and_associations/pombase_phenotype_annotation.phaf.tsv"
    fypo_slim_ids_and_names_table: List[Path] = field(default_factory=lambda: [POMBASE_DATA_PATH / "ontologies_and_associations/fypo_slim_ids_and_names.tsv"])

    def validate_paths(self) -> None:
        """Validate that all required files exist."""
        required_files = [
            self.fypo_obo,
            self.fypo_gaf,
            *self.fypo_slim_ids_and_names_table
        ]
        
        for file_path in required_files:
            if not file_path.exists():
                raise FileNotFoundError(f"FYPO data file not found: {file_path}")
            if not file_path.is_file():
                raise ValueError(f"Path is not a file: {file_path}")

@dataclass
class MONDODataConfig:
    """Configuration for Mondo Disease Ontology data file paths."""
    
    mondo_obo: Path = POMBASE_DATA_PATH / "ontologies_and_associations/mondo-simple.obo"
    mondo_gaf: Path = POMBASE_DATA_PATH / "ontologies_and_associations/human_disease_association.tsv"
    mondo_slim_ids_and_names_table: List[Path] = field(default_factory=lambda: [POMBASE_DATA_PATH / "ontologies_and_associations/pombe_mondo_disease_slim_terms.tsv"])
    
    def validate_paths(self) -> None:
        """Validate that all required files exist."""
        required_files = [
            self.mondo_obo,
            self.mondo_gaf,
            *self.mondo_slim_ids_and_names_table
        ]
        
        for file_path in required_files:
            if not file_path.exists():
                raise FileNotFoundError(f"MONDO data file not found: {file_path}")
            if not file_path.is_file():
                raise ValueError(f"Path is not a file: {file_path}")
    

# ================================= Main Configuration Class =================================

@dataclass
class DataConfig:
    """Main configuration class for all data files."""
    
    gene_metadata: GeneMetadataConfig
    insertion_level: InsertionLevelConfig
    gene_level: GeneLevelConfig
    gene_ontology_data: GeneOntologyDataConfig
    gene_phenotype_data: FYPODataConfig
    mondo_disease_ontology_data: MONDODataConfig

    def validate_all_paths(self) -> None:
        """Validate all file paths across all configurations."""
        self.gene_metadata.validate_paths()
        self.insertion_level.validate_paths()
        self.gene_level.validate_paths()
        self.gene_ontology_data.validate_paths()
        self.gene_phenotype_data.validate_paths()
        self.mondo_disease_ontology_data.validate_paths()

    @classmethod
    def create_default(cls) -> "DataConfig":
        """Create a default configuration with standard file paths."""
        return cls(
            gene_metadata=GeneMetadataConfig(),
            insertion_level=InsertionLevelConfig(),
            gene_level=GeneLevelConfig(),
            gene_ontology_data=GeneOntologyDataConfig(),
            gene_phenotype_data=FYPODataConfig(),
            mondo_disease_ontology_data=MONDODataConfig()
        )


# ================================= Utility Functions =================================

def get_default_config() -> DataConfig:
    """Get the default data configuration."""
    return DataConfig.create_default()

def get_long_timecourse_config() -> DataConfig:
    """Get the long timecourse data configuration."""
    config = DataConfig.create_default()
    config.insertion_level.insertion_annotations = LONG_TIMECOURSE_DATA_PATH / "insertion_level/annotations.tsv.gz"
    config.insertion_level.insertion_level_baseMean = LONG_TIMECOURSE_DATA_PATH / "insertion_level/baseMean.tsv"
    config.insertion_level.insertion_level_LFCs = LONG_TIMECOURSE_DATA_PATH / "insertion_level/LFC.tsv"
    config.insertion_level.insertion_level_weights = LONG_TIMECOURSE_DATA_PATH / "insertion_level/transformed_weights.tsv"
    config.insertion_level.insertion_level_fitting_LFCs = LONG_TIMECOURSE_DATA_PATH / "insertion_level/fitting_LFCs.tsv"
    config.insertion_level.insertion_level_fitting_results = LONG_TIMECOURSE_DATA_PATH / "insertion_level/fitting_results.tsv"
    config.insertion_level.imputation_statistics = None
    config.gene_level.gene_level_LFCs = LONG_TIMECOURSE_DATA_PATH / "gene_level/LFC.tsv"
    config.gene_level.gene_level_fitting_LFCs = LONG_TIMECOURSE_DATA_PATH / "gene_level/fitting_LFCs.tsv"
    config.gene_level.gene_level_fitting_results = LONG_TIMECOURSE_DATA_PATH / "gene_level/fitting_results.tsv"

    return config

def get_haploid_config() -> DataConfig:
    """Get the haploid data configuration."""
    config = DataConfig.create_default()
    config.insertion_level.insertion_annotations = HAPLOID_DATA_PATH / "insertion_level/annotations.tsv.gz"
    config.insertion_level.insertion_level_baseMean = HAPLOID_DATA_PATH / "insertion_level/baseMean.tsv"
    config.insertion_level.insertion_level_LFCs = HAPLOID_DATA_PATH / "insertion_level/LFC.tsv"
    config.insertion_level.insertion_level_weights = HAPLOID_DATA_PATH / "insertion_level/transformed_weights.tsv"
    config.insertion_level.insertion_level_fitting_LFCs = HAPLOID_DATA_PATH / "insertion_level/fitting_LFCs.tsv"
    config.insertion_level.insertion_level_fitting_results = HAPLOID_DATA_PATH / "insertion_level/fitting_results.tsv"
    config.insertion_level.imputation_statistics = None
    config.gene_level.gene_level_LFCs = HAPLOID_DATA_PATH / "gene_level/LFC.tsv"
    config.gene_level.gene_level_fitting_LFCs = HAPLOID_DATA_PATH / "gene_level/fitting_LFCs.tsv"
    config.gene_level.gene_level_fitting_results = HAPLOID_DATA_PATH / "gene_level/fitting_results.tsv"

    return config

def get_custom_config(**kwargs) -> DataConfig:
    """Get a custom data configuration."""
    config = DataConfig.create_default()
    for key, value in kwargs.items():
        setattr(config, key, value)
    return config

def validate_config(config: DataConfig) -> None:
    """Validate a data configuration."""
    config.validate_all_paths()


# ================================= Configuration Instance =================================

# Default configuration instance
DEFAULT_CONFIG = get_default_config()