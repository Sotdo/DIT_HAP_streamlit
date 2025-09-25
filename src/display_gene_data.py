"""
Display gene data.
"""

# ================================= Imports =================================
import numpy as np
import pandas as pd
import altair as alt
import streamlit as st

# ================================= Constants =================================
GENE_BODY_Y = -0.5

# ================================= Functions =================================
def fitting_function(x: np.ndarray, A: float, um: float, lam: float) -> np.ndarray:
    """Calculate Gompertz function values with numerical stability."""
    if A == 0:
        return np.zeros_like(x)
    exponent = np.clip((um * np.e / A) * (lam - x) + 1, -700, 700)
    return A * np.exp(-np.exp(exponent))

@st.cache_data
def display_basic_information(gene: str, gene_info: pd.DataFrame) -> None:
    
    name = gene_info["gene_name"].values[0]
    chromosome = gene_info["chromosome_id"].values[0]
    gene_product = gene_info["gene_product"].values[0]
    uniprot_id = gene_info["uniprot_id"].values[0]
    gene_type = gene_info["gene_type"].values[0]
    synonyms = gene_info["synonyms"].values[0]
    essentiality = gene_info["Gene dispensability. This study"].values[0]
    deletion_mutant_phenotype_description = gene_info["Deletion mutant phenotype description"].values[0]
    phenotypic_classification_used_for_analysis = gene_info["Phenotypic classification used for analysis"].values[0]
    category = gene_info["Category"].values[0]

    pombase_link = f"https://www.pombase.org/gene/{gene}"

    if gene == name:
        st.header(f"Gene: [{gene}]({pombase_link})")
    else:
        st.header(f"Gene: [{name} / {gene}]({pombase_link})")
    st.divider()
    info_dict = {
        "Chromosome": chromosome,
        "Feature type": gene_type,
        "Synonyms": synonyms,
        "Product": gene_product,
        "Uniprot ID": uniprot_id,
        "Essentiality (Hayles 2013)": essentiality,
        "Deletion phenotype": deletion_mutant_phenotype_description,
        "Phenotypic classification used for analysis": phenotypic_classification_used_for_analysis,
        "Category": category
    }

    html_table = "<table>"
    for info_name, info_value in info_dict.items():
        html_table += f"<tr><td style='vertical-align: top; width: 320px;'><b>{info_name}</b></td><td style='vertical-align: top;'>{info_value}</td></tr>"
    html_table += "</table>"
    st.html(html_table)
    st.divider()

def display_gene_body(gene_body: pd.DataFrame) -> tuple[int, alt.Chart]:
    gene_body_copy = gene_body.copy()
    gene_body_copy["um"] = GENE_BODY_Y
    gene_length =  gene_body_copy["ParentalRegion_length"].values[0]

    gene_body_plot = alt.Chart(gene_body_copy).mark_rule().encode(
        x=alt.X("Start_Fraction_to_start_codon:Q", title="Fraction to start codon"),
        x2=alt.X2("End_Fraction_to_start_codon:Q", title=""),
        y=alt.Y("um:Q", title=""),
        color=alt.ColorValue("blue"),
        size=alt.Size("Feature:N", title="Feature", scale=alt.Scale(range=[6, 1], domain=["CDS", "intron"]), legend=None),
    ) + alt.Chart(pd.DataFrame({"x": [1], "y": [GENE_BODY_Y]})).mark_point(shape="triangle", filled=True, fillOpacity=1).encode(
        x=alt.X("x:Q"),
        y=alt.Y("y:Q"),
        angle=alt.AngleValue(90),
        size=alt.SizeValue(150),
        color=alt.ColorValue("blue")
    )

    return gene_length, gene_body_plot

def display_insertion_level_data(gene_length: int, insertion_level_anno_and_results: pd.DataFrame, insertion_level_data: pd.DataFrame, timepoints: dict) -> alt.Chart:

    point_selector = alt.selection_point(fields=["Chr", "Coordinate", "Strand", "Target"], empty=True)

    insertion_level_data_plot1 = alt.Chart(insertion_level_anno_and_results.reset_index()).mark_circle(opacity=0.6).encode(
        x=alt.X("Distance_to_start_codon:Q", title="Distance to start codon", scale=alt.Scale(domain=(0, gene_length))),
        y=alt.Y("um:Q", title="DR", scale=alt.Scale(domain=(-0.5, 1.5))),
        color=alt.condition(point_selector, "Insertion_direction:N", alt.value("lightgray"), legend=alt.Legend(orient="right", title="Insertion direction")),
        size=alt.Size("baseMean:Q", title="baseMean", scale=alt.Scale(type='pow')),
        tooltip=insertion_level_anno_and_results.columns.tolist()
    ).add_params(point_selector)

    insertion_level_data_plot2 = alt.Chart(insertion_level_data.reset_index()).mark_circle(opacity=0.6).encode(
        x=alt.X("Generations:Q", title="Generations", scale=alt.Scale(domain=(0, int(max(timepoints.values()))+1))),
        y=alt.Y("LFC:Q", title="LFC", scale=alt.Scale(domain=(-3, 8))),
        color=alt.condition(point_selector, "num_of_imputed_insertions:N", alt.value("lightgray"), legend=alt.Legend(orient="right", title="Imputation level")),
        size=alt.Size("weights:Q", title="-log10(padj)"),
    ).add_params(point_selector).transform_filter(point_selector)

    return insertion_level_data_plot1, insertion_level_data_plot2

@st.cache_data
def display_gene_level_data(gene_level_fitting_results_in_current_gene: pd.DataFrame, gene_level_data: pd.DataFrame, timepoints: dict) -> alt.Chart:

    fitting_result = gene_level_fitting_results_in_current_gene.copy()
    fitting_result["x"] = 0
    fitting_result["x2"] = 1
    
    gene_level_DR_line = alt.Chart(fitting_result).mark_line(size=4, color="red", opacity=0.3).encode(
        x=alt.X("x:Q", title=""),
        x2=alt.X2("x2:Q", title=""),
        y=alt.Y("um:Q", title=""),
        tooltip=fitting_result.columns.tolist()
    )

    gene_level_data_plot = alt.Chart(gene_level_data.reset_index()).mark_circle(opacity=0.6).encode(
        x=alt.X("Generations:Q", title="", scale=alt.Scale(domain=(0, int(max(timepoints.values()))+1))),
        y=alt.Y("LFC:Q", title="", scale=alt.Scale(domain=(-3, 8))),
        size=alt.value(100),
        color=alt.value("lightgray"),
        tooltip=gene_level_data.columns.tolist()
    )

    A = fitting_result["A"].values[0]
    um = fitting_result["um"].values[0]
    lam = fitting_result["lam"].values[0]

    DL_line = pd.DataFrame(
        {
            "x": [0],
            "x2": [lam],
            "y": [0]
        }
    )

    DL_line_plot = alt.Chart(DL_line).mark_rule(color="green").encode(
        x=alt.X("x:Q", title=""),
        x2=alt.X2("x2:Q", title=""),
        y=alt.Y("y:Q", title="")
    )

    DR_line_x = np.linspace(lam, min(lam + 7/abs(um), int(max(timepoints.values()))+1), 100)
    DR_line_y = (DR_line_x-lam)*um
    DR_line = pd.DataFrame(
        {
            "x": DR_line_x, 
            "y": DR_line_y
        }
    )
    DR_line_plot = alt.Chart(DR_line).mark_line(color="blue").encode(
        x=alt.X("x:Q", title=""),
        y=alt.Y("y:Q", title=""),
    )

    fitting_curve = pd.DataFrame(
        {
            "x": np.linspace(0, int(max(timepoints.values()))+1, 100), 
            "y": fitting_function(np.linspace(0, int(max(timepoints.values()))+1, 100), A, um, lam)
        }
    )

    fitting_curve_plot = alt.Chart(fitting_curve).mark_line(color="red").encode(
        x=alt.X("x:Q", title=""),
        y=alt.Y("y:Q", title=""),
    )

    return gene_level_DR_line, gene_level_data_plot, DL_line_plot, DR_line_plot, fitting_curve_plot


def combine_plots(
    gene_body_plot: alt.Chart,
    insertion_level_data_plot1: alt.Chart,
    insertion_level_data_plot2: alt.Chart,
    gene_level_DR_line: alt.Chart,
    gene_level_data_plot: alt.Chart,
    DL_line_plot: alt.Chart,
    DR_line_plot: alt.Chart,
    fitting_curve_plot: alt.Chart,
) -> alt.Chart:

    plot1 = alt.layer(
        alt.layer(
            gene_body_plot,
            gene_level_DR_line
        ).resolve_scale(color='independent', size="independent"),
        insertion_level_data_plot1
    ).resolve_scale(x='independent', color='independent', size='independent')

    plot2 = alt.layer(
        insertion_level_data_plot2,
        gene_level_data_plot,
        DL_line_plot,
        DR_line_plot,
        fitting_curve_plot
    )

    combined_plot = alt.vconcat(
        plot1.properties(height=200),
        plot2.properties(height=200)
    ).resolve_scale(x='independent', y='independent', size='independent', color='independent')

    return combined_plot

def display_gene_level_metrics(container: st.container, gene_level_fitting_results_in_current_gene: pd.DataFrame):

    with container:
        st.metric(label="DR", value=round(gene_level_fitting_results_in_current_gene["um"].values[0], 3))
        st.metric(label="DL", value=round(gene_level_fitting_results_in_current_gene["lam"].values[0], 3))
        st.metric(label="R^2", value=round(gene_level_fitting_results_in_current_gene["R2"].values[0], 3))
        st.metric(label="RMSE", value=round(gene_level_fitting_results_in_current_gene["RMSE"].values[0], 3))
        st.metric(label="Normalized RMSE", value=round(gene_level_fitting_results_in_current_gene["normalized_RMSE"].values[0], 3))