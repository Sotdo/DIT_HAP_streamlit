"""
Enrichment analysis page using the new modular data loading architecture.

This page provides gene ontology and phenotype enrichment analysis functionality.
"""

# ================================= Imports =================================
import streamlit as st
import sys
import pandas as pd
sys.path.append("../src")
from src.preparation import sidebar_gene_input
from src.data_config import get_default_config
from src.data_manager import load_gene_metadata, load_gene_level_stats, load_gene_ontology_data, load_gene_phenotype_data, load_disease_ontology_data
from src.enrichment_functions import ontology_enrichment_pipeline, stringdb_enrichment, display_enrichment_results

# ================================= Constants =================================
P_VALUE_THRESHOLD = 0.05

# ================================= Functions =================================

def load_enrichment_data():
    """Load only the data needed for enrichment analysis."""
    
    # Get default configuration
    config = get_default_config()
    
    # Validate configuration
    config.validate_all_paths()
    
    # Load only the required data categories
    with st.spinner("Loading gene metadata...", show_time=True):
        gene_metadata = load_gene_metadata(config.gene_metadata)

    with st.spinner("Loading gene level statistics...", show_time=True):
        gene_level = load_gene_level_stats(config.gene_level)
    
    with st.spinner("Loading ontology data...", show_time=True):
        gene_ontology_data = load_gene_ontology_data(config.gene_ontology_data)
    
    with st.spinner("Loading phenotype data...", show_time=True):
        gene_phenotype_data = load_gene_phenotype_data(config.gene_phenotype_data)
    
    with st.spinner("Loading disease ontology data...", show_time=True):
        disease_ontology_data = load_disease_ontology_data(config.mondo_disease_ontology_data)
    
    return gene_metadata, gene_level, gene_ontology_data, gene_phenotype_data, disease_ontology_data

def display_results(res: pd.DataFrame, res_slim: pd.DataFrame = None):
    """Display enrichment results."""

    
    if res.empty:
        st.warning("No enrichment results found")
    else:
        st.success("Enrichment results found")
        st.altair_chart(display_enrichment_results(res), use_container_width=True)

    if res_slim is not None:
        if res_slim.empty:
            st.warning("No enrichment results found (slim)")
        else:
            st.success("Enrichment results found (slim)")
            st.altair_chart(display_enrichment_results(res_slim), use_container_width=True)

def main():
    """Main entry point for the enrichment analysis page."""
    
    # Load required data
    gene_metadata, gene_level, gene_ontology_data, gene_phenotype_data, disease_ontology_data = load_enrichment_data()
    
    # Get gene input from sidebar
    covered_gene_sysIDs, submit_button = sidebar_gene_input(
        gene_metadata.gene_info_with_essentiality, 
        gene_level.gene_level_LFCs
    )
    
    if submit_button and covered_gene_sysIDs:

        bg_genes = gene_level.gene_level_LFCs.index.tolist()
        
        ontology_tab = st.tabs(["GO enrichment", "FYPO enrichment", "Mondo enrichment", "STRING enrichment"])
        with ontology_tab[0]:
            with st.spinner("Performing GO enrichment analysis..."):

                load_kwargs = {
                    "relationships": {"is_a", "part_of"},
                    "propagate_counts": True,
                    "load_obsolete": False
                }
                
                enrichment_kwargs = {
                    "alpha": P_VALUE_THRESHOLD,
                    "methods": ["fdr_bh"],
                    "propagate_counts": True,
                    "relationships": {"is_a", "part_of"},
                    "prt": None,
                }

                format_kwargs = {
                    "itemid2name": gene_metadata.id2name
                }

                res, res_slim = ontology_enrichment_pipeline(gene_ontology_data, covered_gene_sysIDs, bg_genes, load_kwargs=load_kwargs, enrichment_kwargs=enrichment_kwargs, format_kwargs=format_kwargs)
                display_results(res, res_slim)
        with ontology_tab[1]:
            with st.spinner("Performing FYPO enrichment analysis..."):

                load_kwargs = {
                    "propagate_counts": True,
                    "load_obsolete": False
                }

                enrichment_kwargs = {
                    "alpha": P_VALUE_THRESHOLD,
                    "methods": ["fdr_bh"],
                    "propagate_counts": True,
                    "prt": None,
                }

                format_kwargs = {
                    "itemid2name": gene_metadata.id2name
                }

                res, res_slim = ontology_enrichment_pipeline(gene_phenotype_data, covered_gene_sysIDs, bg_genes, load_kwargs=load_kwargs, enrichment_kwargs=enrichment_kwargs, format_kwargs=format_kwargs)
                display_results(res, res_slim)
        with ontology_tab[2]:
            with st.spinner("Performing Mondo enrichment analysis..."):

                load_kwargs = {
                    "propagate_counts": True,
                    "load_obsolete": False
                }

                enrichment_kwargs = {
                    "alpha": P_VALUE_THRESHOLD,
                    "methods": ["fdr_bh"],
                    "propagate_counts": True,
                    "prt": None,
                }

                format_kwargs = {
                    "itemid2name": gene_metadata.id2name
                }
                res, res_slim = ontology_enrichment_pipeline(disease_ontology_data, covered_gene_sysIDs, bg_genes, load_kwargs=load_kwargs, enrichment_kwargs=enrichment_kwargs, format_kwargs=format_kwargs)
                display_results(res, res_slim)
        with ontology_tab[3]:
            with st.spinner("Performing STRING enrichment analysis..."):
                res = stringdb_enrichment(covered_gene_sysIDs, bg_genes)
                display_results(res)
    else:
        st.warning("Please select genes for enrichment analysis")


if __name__ == "__main__":
    main()