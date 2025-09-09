"""
Extract DIT-HAP data for a given gene.
"""

# ================================= Imports =================================
from pathlib import Path
from typing import List
from dataclasses import dataclass
import numpy as np
import pandas as pd
import streamlit as st
from .data_manager import GeneMetadataData, InsertionLevelData, GeneLevelData


# ================================= Functions =================================

def get_gene_info(gene: str, gene_info_with_essentiality: pd.DataFrame) -> pd.DataFrame:
    """Get the gene info for a given gene."""
    gene_info = gene_info_with_essentiality.query(
        "gene_systematic_id == @gene"
    )
    return gene_info

def get_gene_body(gene: str, genome_intervals: pd.DataFrame) -> pd.DataFrame:
    """Get the gene body for a given gene."""
    gene_body = genome_intervals.query(
        "`Systematic ID` == @gene"
    )
    return gene_body


def get_insertion_level_data(
    gene: str, 
    insertion_level_data: InsertionLevelData,
    timepoints: dict
) -> pd.DataFrame:
    """Get the insertions in a given gene."""

    tps = list(timepoints.keys())

    insertions_in_current_gene = insertion_level_data.insertion_annotations_with_baseMean.query(
        "(`Systematic ID` == @gene) and (Distance_to_stop_codon > 4)"
    ).index.intersection(insertion_level_data.insertion_level_LFCs.index)

    annotation_in_current_gene = insertion_level_data.insertion_annotations_with_baseMean.loc[insertions_in_current_gene]
    insertion_LFCs_in_current_gene = insertion_level_data.insertion_level_LFCs.loc[insertions_in_current_gene][tps].rename_axis("Timepoint", axis=1).stack().rename("LFC")
    insertion_weights_in_current_gene = insertion_level_data.insertion_level_weights.loc[insertions_in_current_gene][tps].rename_axis("Timepoint", axis=1).stack().rename("weights")
    insertion_fitting_LFCs_in_current_gene = insertion_level_data.insertion_level_fitting_LFCs.loc[insertions_in_current_gene][tps].rename_axis("Timepoint", axis=1).stack().rename("fitting_LFC")
    insertion_fitting_results_in_current_gene = insertion_level_data.insertion_level_fitting_results.loc[insertions_in_current_gene]
    insertion_level_data_df = pd.concat([insertion_LFCs_in_current_gene, insertion_weights_in_current_gene, insertion_fitting_LFCs_in_current_gene], axis=1)
    insertion_level_data_df["Generations"] = [ timepoints[tp] for tp in insertion_level_data_df.index.get_level_values(-1)]
    insertion_level_data_df = insertion_level_data_df.reset_index().set_index(["Chr", "Coordinate", "Strand", "Target"])

    insertion_level_anno_and_results = insertion_fitting_results_in_current_gene.merge(
        annotation_in_current_gene,
        how="left",
        left_index=True,
        right_index=True,
    )

    insertion_level_data_df = insertion_level_data_df.merge(
        insertion_level_anno_and_results,
        how="left",
        left_index=True,
        right_index=True,
    )

    return insertion_level_anno_and_results, insertion_level_data_df

def get_gene_level_data(
    gene: str, 
    gene_level_data: GeneLevelData,
    timepoints: dict
) -> pd.DataFrame:
    """Get the gene level data for a given gene."""

    tps = list(timepoints.keys())
    gene_level_LFCs_in_current_gene = gene_level_data.gene_level_LFCs.loc[[gene], tps].rename_axis("Timepoint", axis=1).stack().to_frame("LFC")
    gene_level_fitting_LFCs_in_current_gene = gene_level_data.gene_level_fitting_LFCs.loc[[gene], tps].rename_axis("Timepoint", axis=1).stack().to_frame("fitting_LFC")
    gene_level_fitting_results_in_current_gene = gene_level_data.gene_level_fitting_results.loc[[gene]]
    gene_level_data_df = pd.concat([gene_level_LFCs_in_current_gene, gene_level_fitting_LFCs_in_current_gene], axis=1)
    gene_level_data_df["Generations"] = [ timepoints[tp] for tp in gene_level_data_df.index.get_level_values(-1)]
    gene_level_data_df = gene_level_data_df.merge(
        gene_level_fitting_results_in_current_gene,
        how="left",
        left_index=True,
        right_index=True,
    )

    return gene_level_fitting_results_in_current_gene, gene_level_data_df
