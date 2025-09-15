"""
Feature space page for DIT-HAP Streamlit application.

This page provides a visualization of the query genes in the feature space.
"""

# ================================= Imports =================================
import streamlit as st
import sys
import pandas as pd
import altair as alt
sys.path.append("../src")
from src.preparation import sidebar_gene_input
from src.data_config import get_default_config
from src.data_manager import load_gene_metadata, load_gene_level_stats, load_gene_ontology_data, load_gene_phenotype_data, load_disease_ontology_data

# ================================= Functions =================================

def load_data():
    """Load only the data needed for feature space analysis."""
    
    # Get default configuration
    config = get_default_config()
    
    # Validate configuration
    config.validate_all_paths()
    
    # Load only the required data categories
    with st.spinner("Loading gene metadata...", show_time=True):
        gene_metadata = load_gene_metadata(config.gene_metadata)

    with st.spinner("Loading gene level statistics...", show_time=True):
        gene_level = load_gene_level_stats(config.gene_level)
    
    return gene_metadata, gene_level

def display_feature_space(query_genes: list[str], gene_level: GeneLevelData) -> alt.Chart:
    """Display the feature space for the query genes."""

    all_gene_feature_space = alt.Chart(gene_level.gene_level_LFCs).mark_circle(opacity=0.6).encode(
        x=alt.X("um:Q", title="Depletion rate"),
        y=alt.Y("lam:Q", title="Depletion lag"),
        color=alt.value("lightgray"),
        tooltip=gene_level.gene_level_LFCs.columns.tolist()
    )

    query_gene_feature_space = alt.Chart(gene_level.gene_level_LFCs.loc[query_genes]).mark_circle(opacity=0.6).encode(
        x=alt.X("um:Q", title="Depletion rate"),
        y=alt.Y("lam:Q", title="Depletion lag"),
        color=alt.value("red"),
        tooltip=gene_level.gene_level_LFCs.columns.tolist()
    )

    return all_gene_feature_space + query_gene_feature_space

def main():
    """Main entry point for the feature space page."""
    
    # Load required data
    gene_metadata, gene_level = load_data()
    
    # Get gene input from sidebar
    covered_gene_sysIDs, submit_button = sidebar_gene_input(
        gene_metadata.gene_info_with_essentiality, 
        gene_level.gene_level_LFCs
    )
    
    if submit_button and covered_gene_sysIDs:

        bg_genes = gene_level.gene_level_LFCs.index.tolist()

        alt_chart = display_feature_space(covered_gene_sysIDs, gene_level)
        st.altair_chart(alt_chart, use_container_width=True)

if __name__ == "__main__":
    main()