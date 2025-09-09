"""
Data manager module for DIT-HAP Streamlit application.

This module provides modular data loading functions for different data categories,
supporting both individual category loading and full data loading.
"""

# ================================= Imports =================================
from typing import List, Optional
from pathlib import Path
from dataclasses import dataclass
import pandas as pd
import streamlit as st
from .data_config import GeneMetadataConfig, InsertionLevelConfig, GeneLevelConfig, GeneOntologyDataConfig, FYPODataConfig, MONDODataConfig
from .preparation import format_phaf_file, format_mondo_gaf_file

# ================================= Data Classes =================================

@dataclass
class GeneMetadataData:
    """Container for gene metadata data."""
    gene_info_with_essentiality: pd.DataFrame
    genome_intervals: pd.DataFrame
    id2name: dict


@dataclass
class InsertionLevelData:
    """Container for insertion level statistics data."""
    insertion_annotations_with_baseMean: pd.DataFrame
    insertion_level_LFCs: pd.DataFrame
    insertion_level_weights: pd.DataFrame
    insertion_level_fitting_LFCs: pd.DataFrame
    insertion_level_fitting_results: pd.DataFrame


@dataclass
class GeneLevelData:
    """Container for gene level statistics data."""
    gene_level_LFCs: pd.DataFrame
    gene_level_fitting_LFCs: pd.DataFrame
    gene_level_fitting_results: pd.DataFrame


@dataclass
class GeneOntologyData:
    """Container for ontology data."""
    ontology_obo: Path
    ontology_association_gaf: Path
    slim_terms: pd.DataFrame

@dataclass
class FYPOData:
    """Container for fission yeast phenotype ontology data."""
    ontology_obo: Path
    ontology_association_gaf: Path
    slim_terms: pd.DataFrame

@dataclass
class MONDOData:
    """Container for Mondo Disease Ontology data."""
    ontology_obo: Path
    ontology_association_gaf: Path
    slim_terms: pd.DataFrame

@dataclass
class DataManager:
    """Complete data manager containing all data categories."""
    gene_metadata: GeneMetadataData
    insertion_level: InsertionLevelData
    gene_level: GeneLevelData
    gene_ontology_data: GeneOntologyData
    gene_phenotype_data: FYPOData
    mondo_disease_ontology_data: MONDOData


# ================================= Utility Functions =================================

def read_file(file: Path, index_col: Optional[List[int]] = None, header: Optional[List[int] | str] = "infer", **kwargs) -> pd.DataFrame:
    """Read a file into a pandas DataFrame based on file extension."""
    if "tsv" in file.name:
        return pd.read_csv(file, index_col=index_col, header=header, sep="\t", **kwargs)
    elif "bed" in file.name:
        return pd.read_csv(file, index_col=index_col, header=header, sep="\t", **kwargs)
    elif "csv" in file.name:
        return pd.read_csv(file, index_col=index_col, header=header, sep=",", **kwargs)
    elif "xlsx" in file.name:
        if header == "infer":
            header = 0
        return pd.read_excel(file, index_col=index_col, header=header, **kwargs)
    else:
        raise ValueError(f"Unsupported file type: {file.name}")


# ================================= Gene Metadata Loading =================================

@st.cache_data
def load_gene_metadata(config: GeneMetadataConfig) -> GeneMetadataData:
    """Load gene metadata data from configuration."""
    
    # Read individual files
    gene_IDs_names_products = read_file(config.gene_IDs_names_products)
    deletion_library_essentiality = read_file(config.deletion_library_essentiality)
    genome_intervals = read_file(config.genome_intervals)

    # Fill NA in gene_name with gene_systematic_id
    gene_IDs_names_products["gene_name"] = gene_IDs_names_products["gene_name"].fillna(gene_IDs_names_products["gene_systematic_id"])
    id2name = dict(zip(list(gene_IDs_names_products["gene_systematic_id"]), list(gene_IDs_names_products["gene_name"])))

    # Merge gene_IDs_names_products and deletion_library_essentiality
    gene_info_with_essentiality = gene_IDs_names_products.merge(
        deletion_library_essentiality[["Updated_Systematic_ID", "Gene dispensability. This study", "Deletion mutant phenotype description", "Phenotypic classification used for analysis", "Category"]],
        how="left",
        left_on="gene_systematic_id",
        right_on="Updated_Systematic_ID"
    ).drop(
        columns=["Updated_Systematic_ID"]
    )

    # Calculate genome interval fractions
    genome_intervals["Start_Fraction_to_region_start"] = (genome_intervals["Start"] - genome_intervals["ParentalRegion_start"]) / genome_intervals["ParentalRegion_length"]
    genome_intervals["End_Fraction_to_region_start"] = (genome_intervals["End"] - genome_intervals["ParentalRegion_start"]) / genome_intervals["ParentalRegion_length"]
    genome_intervals["Start_Fraction_to_region_end"] = (genome_intervals["ParentalRegion_end"] - genome_intervals["Start"]) / genome_intervals["ParentalRegion_length"]
    genome_intervals["End_Fraction_to_region_end"] = (genome_intervals["ParentalRegion_end"] - genome_intervals["End"]) / genome_intervals["ParentalRegion_length"]

    def fraction_to_start_codon(row):
        strand = row["Strand"]
        if strand == "+":
            return row["Start_Fraction_to_region_start"], row["End_Fraction_to_region_start"]
        else:
            return row["End_Fraction_to_region_end"], row["Start_Fraction_to_region_end"]

    genome_intervals[["Start_Fraction_to_start_codon", "End_Fraction_to_start_codon"]] = genome_intervals.apply(fraction_to_start_codon, axis=1, result_type="expand")
    
    return GeneMetadataData(
        gene_info_with_essentiality=gene_info_with_essentiality,
        genome_intervals=genome_intervals,
        id2name=id2name
    )


# ================================= Insertion Level Loading =================================

@st.cache_data
def load_insertion_level_stats(config: InsertionLevelConfig) -> InsertionLevelData:
    """Load insertion level statistics data from configuration."""
    
    index_col = [0, 1, 2, 3]  # chr, coordinate, strand, target
    
    # Read all insertion level files
    insertion_annotations = read_file(config.insertion_annotations, index_col=index_col)
    insertion_level_baseMean = read_file(config.insertion_level_baseMean, index_col=index_col).iloc[:, 0].rename("baseMean")
    insertion_level_LFCs = read_file(config.insertion_level_LFCs, index_col=index_col)
    insertion_level_weights = read_file(config.insertion_level_weights, index_col=index_col)
    insertion_level_fitting_LFCs = read_file(config.insertion_level_fitting_LFCs, index_col=index_col)
    insertion_level_fitting_results = read_file(config.insertion_level_fitting_results, index_col=index_col)

    # Merge insertion_annotations and insertion_level_baseMean
    insertion_annotations_with_baseMean = insertion_annotations.merge(
        insertion_level_baseMean,
        how="left",
        left_index=True,
        right_index=True,
    )

    # Merge with imputation statistics if available
    if config.imputation_statistics and config.imputation_statistics.exists():
        imputation_statistics = read_file(config.imputation_statistics, index_col=index_col)
        insertion_annotations_with_baseMean = insertion_annotations_with_baseMean.merge(
            imputation_statistics,
            how="left",
            left_index=True,
            right_index=True,
        )

    return InsertionLevelData(
        insertion_annotations_with_baseMean=insertion_annotations_with_baseMean,
        insertion_level_LFCs=insertion_level_LFCs,
        insertion_level_weights=insertion_level_weights,
        insertion_level_fitting_LFCs=insertion_level_fitting_LFCs,
        insertion_level_fitting_results=insertion_level_fitting_results
    )


# ================================= Gene Level Loading =================================

@st.cache_data
def load_gene_level_stats(config: GeneLevelConfig) -> GeneLevelData:
    """Load gene level statistics data from configuration."""
    
    index_col = 0
    gene_level_LFCs = read_file(config.gene_level_LFCs, index_col=index_col)
    gene_level_fitting_LFCs = read_file(config.gene_level_fitting_LFCs, index_col=index_col)
    gene_level_fitting_results = read_file(config.gene_level_fitting_results, index_col=index_col)

    return GeneLevelData(
        gene_level_LFCs=gene_level_LFCs,
        gene_level_fitting_LFCs=gene_level_fitting_LFCs,
        gene_level_fitting_results=gene_level_fitting_results
    )


# ================================= Ontology Data Loading =================================

@st.cache_data
def load_gene_ontology_data(config: GeneOntologyDataConfig) -> GeneOntologyData:
    """Load ontology data configuration (files are loaded by enrichment functions)."""

    goslim_terms = pd.concat([read_file(path, header=None, names=["Term", "Description"]) for path in config.goslim_terms_table])
    
    return GeneOntologyData(
        ontology_obo=config.gene_ontology_obo,
        ontology_association_gaf=config.gene_ontology_association_gaf,
        slim_terms=goslim_terms
    )

@st.cache_data
def load_gene_phenotype_data(config: FYPODataConfig) -> FYPOData:
    """Load gene phenotype data from configuration."""
    fypo_slim_ids_and_names = pd.concat([read_file(path, header=None, names=["Term", "Description"]) for path in config.fypo_slim_ids_and_names_table])
    return FYPOData(
        ontology_obo=config.fypo_obo,
        ontology_association_gaf=format_phaf_file(config.fypo_obo, config.fypo_gaf),
        slim_terms=fypo_slim_ids_and_names
    )

@st.cache_data
def load_disease_ontology_data(config: MONDODataConfig) -> MONDOData:
    """Load disease ontology data from configuration."""

    mondo_slim_ids_and_names = pd.concat([read_file(path, header=None, names=["Term", "Description"]) for path in config.mondo_slim_ids_and_names_table])
    
    return MONDOData(
        ontology_obo=config.mondo_obo,
        ontology_association_gaf=format_mondo_gaf_file(config.mondo_obo, config.mondo_gaf),
        slim_terms=mondo_slim_ids_and_names
    )

# ================================= Complete Data Loading =================================

@st.cache_data
def load_all_data(
    gene_config: GeneMetadataConfig,
    insertion_config: InsertionLevelConfig, 
    gene_level_config: GeneLevelConfig, 
    gene_ontology_data_config: GeneOntologyDataConfig, 
    gene_phenotype_data_config: FYPODataConfig, 
    mondo_disease_ontology_data_config: MONDODataConfig
) -> DataManager:
    """Load all data categories into a complete DataManager."""
    
    gene_metadata = load_gene_metadata(gene_config)
    insertion_level = load_insertion_level_stats(insertion_config)
    gene_level = load_gene_level_stats(gene_level_config)
    gene_ontology_data = load_gene_ontology_data(gene_ontology_data_config)
    gene_phenotype_data = load_gene_phenotype_data(gene_phenotype_data_config)
    mondo_disease_ontology_data = load_disease_ontology_data(mondo_disease_ontology_data_config)

    return DataManager(
        gene_metadata=gene_metadata,
        insertion_level=insertion_level,
        gene_level=gene_level,
        gene_ontology_data=gene_ontology_data,
        gene_phenotype_data=gene_phenotype_data,
        mondo_disease_ontology_data=mondo_disease_ontology_data
    )


# ================================= Backward Compatibility =================================

@st.cache_data
def load_data_from_config(config) -> DataManager:
    """Load data using the full configuration (backward compatibility)."""
    
    return load_all_data(
        gene_config=config.gene_metadata,
        insertion_config=config.insertion_level,
        gene_level_config=config.gene_level,
        gene_ontology_data_config=config.gene_ontology_data,
        gene_phenotype_data_config=config.gene_phenotype_data,
        mondo_disease_ontology_data_config=config.mondo_disease_ontology_data
    )