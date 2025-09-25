"""
Updated plot page using the new modular data loading architecture.

This version demonstrates how to use the new simplified data loading
system with modular configuration.
"""

# ================================= Imports =================================
import streamlit as st
import altair as alt
import pandas as pd
import sys
sys.path.append("../src")
from src.preparation import sidebar_gene_input
from src.data_config import get_default_config, get_long_timecourse_config, get_haploid_config, DataConfig
from src.data_manager import load_gene_metadata, load_insertion_level_stats, load_gene_level_stats, GeneMetadataData, InsertionLevelData, GeneLevelData
from src.get_gene_data import get_gene_info, get_gene_body, get_insertion_level_data, get_gene_level_data
from src.display_gene_data import display_basic_information, display_gene_body, display_insertion_level_data, display_gene_level_data, combine_plots, display_gene_level_metrics

# ================================= Constants =================================
P_VALUE_THRESHOLD = 0.05

TIME_POINTS = {
    "YES0": 0,
    "YES1": 3.352,
    "YES2": 6.588,
    "YES3": 10.104,
    "YES4": 13.480
}

TIME_POINTS_LONG_TIMECOURSE = {
    "YES0": 0,
    "YES1": 3.723,
    "YES2": 6.969,
    "YES3": 10.104,
    "YES4": 13.554,
    "YES5": 17.098,
    "YES6": 20.059
}

TIME_POINTS_HAPLOID = {
    "0h": 0,
    "YES0": 0.553,
    "YES1": 2.097,
    "YES2": 5.629,
    "YES3": 8.831,
    "YES4": 12.203,
    "YES5": 15.818,
    "YES6": 19.081
}

# ================================= Functions =================================

def load_depletion_data(config: DataConfig) -> tuple[GeneMetadataData, InsertionLevelData, GeneLevelData]:
    """Load only the data needed for depletion analysis."""
    
    # Validate configuration
    config.validate_all_paths()
    
    # Load only the required data categories
    with st.spinner("Loading gene metadata...", show_time=True):
        gene_metadata = load_gene_metadata(config.gene_metadata)
    
    with st.spinner("Loading insertion level statistics...", show_time=True):
        insertion_level = load_insertion_level_stats(config.insertion_level)
    
    with st.spinner("Loading gene level statistics...", show_time=True):
        gene_level = load_gene_level_stats(config.gene_level)
    
    return gene_metadata, insertion_level, gene_level

def get_gene_result(
    gene: str,
    gene_length: int,
    timepoints: dict,
    insertion_level: InsertionLevelData,
    gene_level: GeneLevelData,
    gene_body_plot: alt.Chart
) -> tuple[pd.DataFrame, alt.Chart, bool]:
    """Get the gene result for a given gene."""
    try:
        insertion_level_anno_and_results, insertion_level_data = get_insertion_level_data(gene, insertion_level, timepoints)
        gene_level_fitting_results_in_current_gene, gene_level_data = get_gene_level_data(gene, gene_level, timepoints)
        
        insertion_level_data_plot1, insertion_level_data_plot2 = display_insertion_level_data(
            gene_length, 
            insertion_level_anno_and_results, 
            insertion_level_data,
            timepoints
        )
        
        gene_level_DR_line, gene_level_data_plot, DL_line_plot, DR_line_plot, fitting_curve_plot = display_gene_level_data(
            gene_level_fitting_results_in_current_gene, 
            gene_level_data,
            timepoints
        )

        combined_plot = combine_plots(
            gene_body_plot, 
            insertion_level_data_plot1, 
            insertion_level_data_plot2, 
            gene_level_DR_line, 
            gene_level_data_plot, 
            DL_line_plot, 
            DR_line_plot, 
            fitting_curve_plot
        )

        return gene_level_fitting_results_in_current_gene, combined_plot, True
    except Exception as e:
        return None, None, False


def main():
    """Main entry point for the depletion data page."""

    # Get default configuration
    config = get_default_config()
    
    # Load required data
    gene_metadata, insertion_level, gene_level = load_depletion_data(config)

    long_timecourse_config = get_long_timecourse_config()
    _, insertion_level_long_timecourse, gene_level_long_timecourse = load_depletion_data(long_timecourse_config)
    
    haploid_config = get_haploid_config()
    _, insertion_level_haploid, gene_level_haploid = load_depletion_data(haploid_config)
    
    # Get gene input from sidebar
    covered_gene_sysIDs, submit_button = sidebar_gene_input(
        gene_metadata.gene_info_with_essentiality, 
        gene_level.gene_level_LFCs
    )

    if submit_button and covered_gene_sysIDs:
        for gene in covered_gene_sysIDs:
            gene_info = get_gene_info(gene, gene_metadata.gene_info_with_essentiality)
            display_basic_information(gene, gene_info)
            
            gene_body = get_gene_body(gene, gene_metadata.genome_intervals)
            gene_length, gene_body_plot = display_gene_body(gene_body) 

            gene_level_fitting_results_in_current_gene, combined_plot, has_data = get_gene_result(
                gene, 
                gene_length, 
                TIME_POINTS,
                insertion_level, 
                gene_level, 
                gene_body_plot
            )

            gene_level_fitting_results_in_current_gene_long_timecourse, combined_plot_long_timecourse, has_data_long_timecourse = get_gene_result(
                gene, 
                gene_length, 
                TIME_POINTS_LONG_TIMECOURSE,
                insertion_level_long_timecourse, 
                gene_level_long_timecourse, 
                gene_body_plot
            )

            gene_level_fitting_results_in_current_gene_haploid, combined_plot_haploid, has_data_haploid = get_gene_result(
                gene, 
                gene_length, 
                TIME_POINTS_HAPLOID,
                insertion_level_haploid, 
                gene_level_haploid, 
                gene_body_plot
            )
            
            col1, col2, col3, col4, col5, col6 = st.columns([2, 8, 2, 8, 2, 8], border=True)
            if has_data:
                with col1:
                    display_gene_level_metrics(col1, gene_level_fitting_results_in_current_gene)
                with col2:
                    st.altair_chart(combined_plot, use_container_width=True, theme=None)
            else:
                with col1:
                    st.warning("No data found")
                with col2:
                    st.warning("No data found")
            if has_data_long_timecourse:
                with col3:
                    display_gene_level_metrics(col3, gene_level_fitting_results_in_current_gene_long_timecourse)
                with col4:
                    st.altair_chart(combined_plot_long_timecourse, use_container_width=True, theme=None)
            else:
                with col3:
                    st.warning("No data found")
                with col4:
                    st.warning("No data found")
            if has_data_haploid:
                with col5:
                    display_gene_level_metrics(col5, gene_level_fitting_results_in_current_gene_haploid)
                with col6:
                    st.altair_chart(combined_plot_haploid, use_container_width=True, theme=None)
            else:
                with col5:
                    st.warning("No data found")
                with col6:
                    st.warning("No data found")

            st.success("Plot generated successfully")
    else:
        st.warning("No genes submitted or no valid genes")


if __name__ == "__main__":
    main()