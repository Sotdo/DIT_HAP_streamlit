"""
Enrichment functions.
"""

# ================================= Imports =================================

import pandas as pd
from .data_manager import GeneOntologyData, FYPOData, MONDOData
import streamlit as st
from pathlib import Path
from typing import Optional, Literal
import time
import requests
from requests.exceptions import ConnectionError, RequestException
from io import StringIO
import altair as alt

# == GOATools imports ==
from goatools.obo_parser import GODag
from goatools.anno.gaf_reader import GafReader
from goatools.goea.go_enrichment_ns import GOEnrichmentStudyNS
from goatools.go_enrichment import GOEnrichmentRecord
from goatools.rpt.goea_nt_xfrm import get_goea_nts_prt

# ================================= Constants =================================
STRING_API_URL = "https://string-db.org/api"
STRING_OUTPUT_FORMAT = Literal["tsv", "tsv-no-header", "json", "xml"]
STRING_METHOD = Literal["get_string_ids", "enrichment"]
STRING_SPECIES_ID = "4896"
STRING_CALLER_IDENTITY = "dit-hap.streamlit.app"
STRING_SEPARATOR = "%0d"
MAX_RETRIES = 5
RETRY_DELAY = 10


# ================================= Helper Functions =================================
def mapslim(term: str, dag: GODag, slim_dag: dict) -> tuple[set[str], set[str]]:
    """Maps a term (accession) to it's slim terms."""

    all_ancestors = set()
    covered_ancestors = set()

    # get all paths for the term in the dag
    paths = dag.paths_to_top(term)
    for path in paths:
        # the next loop needs to run bottom->up, i.e. from the term item to
        # the root, thus we need to reverse the list prior to iteration
        path.reverse()

        got_leaf = False
        for term in path:
            if term.id in slim_dag:
                all_ancestors.add(term.id)
                if got_leaf:
                    covered_ancestors.add(term.id)
                got_leaf = True

    # get the direct ancestors, i.e. those that are not covered by a earlier
    # ancestor of the slim in _any_ path (in bottom->top order)
    direct_ancestors = all_ancestors - covered_ancestors
    return direct_ancestors, all_ancestors


def get_slim_ns2assoc(ns2assoc: dict, dag: GODag, slim_dag: dict) -> dict:
    """Get the slim ns2assoc."""

    term2slim = {}
    for term in dag:
        term2slim[term] = {}
        term2slim[term]["direct_ancestors"], term2slim[term]["all_ancestors"] = mapslim(
            term, dag, slim_dag
        )

    ns2slim_assoc = {"direct_ancestors": {}, "all_ancestors": {}}
    for ns, gene2terms in ns2assoc.items():
        ns2slim_assoc["direct_ancestors"][ns] = {}
        ns2slim_assoc["all_ancestors"][ns] = {}
        for gene, terms in gene2terms.items():
            ns2slim_assoc["direct_ancestors"][ns][gene] = set()
            ns2slim_assoc["all_ancestors"][ns][gene] = set()
            for term in terms:
                ns2slim_assoc["direct_ancestors"][ns][gene].update(
                    term2slim[term]["direct_ancestors"]
                )
                ns2slim_assoc["all_ancestors"][ns][gene].update(
                    term2slim[term]["all_ancestors"]
                )
    return ns2slim_assoc


def create_enrichment_dataframe(
    oea_results_sig: list[GOEnrichmentRecord], **kwargs
) -> pd.DataFrame:
    """Create an enrichment dataframe from the GOEnrichmentRecord objects manually."""
    results_list = []
    for result in oea_results_sig:
        res_dict = {
            "GO": result.GO,
            "NS": result.NS,
            "name": result.name,
            "level": result.goterm.level,
            "depth": result.goterm.depth,
            "p_uncorrected": result.p_uncorrected,
            "p_fdr_bh": result.p_fdr_bh,
            "study_count": result.study_count,
            "study_n": result.study_n,
            "pop_count": result.pop_count,
            "pop_n": result.pop_n,
            "ratio_in_study": "/".join(map(str, result.ratio_in_study)),
            "ratio_in_pop": "/".join(map(str, result.ratio_in_pop)),
            "study_items": result.study_items,
            "pop_items": result.pop_items,
        }

        # get the definition, but handle the case where it is not available
        try:
            res_dict["defn"] = result.goterm.defn
        except Exception as e:
            res_dict["defn"] = ""

        # convert the study and pop items to names, but handle the case where the itemid2name is not available
        if kwargs["itemid2name"] is not None:
            study_items = sorted([
                kwargs["itemid2name"][item] for item in result.study_items
            ])
            pop_items = sorted([
                kwargs["itemid2name"][item] for item in result.pop_items
            ])
            res_dict["study_items"] = ", ".join(study_items)
            res_dict["pop_items"] = ", ".join(pop_items)
        else:
            res_dict["study_items"] = ", ".join(sorted(result.study_items))
            res_dict["pop_items"] = ", ".join(sorted(result.pop_items))

        # append the result to the list
        results_list.append(res_dict)

    oea_results_sig_prt = pd.DataFrame(results_list)
    return oea_results_sig_prt


# ================================= Main Functions =================================
# @st.cache_resource
def load_ontology_data(
    ontology_data: GeneOntologyData | FYPOData | MONDOData, **kwargs
) -> tuple[GODag, GafReader, dict]:
    """Load ontology data from obo file and association file."""
    try:
        dag = GODag(
            str(ontology_data.ontology_obo),
            optional_attrs=["def", "relationship"],
            load_obsolete=False,
        )
    except KeyError:
        dag = GODag(
            str(ontology_data.ontology_obo), optional_attrs=["def"], load_obsolete=False
        )

    objanno = GafReader(str(ontology_data.ontology_association_gaf), godag=dag)

    slim_terms = ontology_data.slim_terms["Term"].to_list()
    slim_dag = {term: dag[term] for term in slim_terms}

    ns2assoc = objanno.get_ns2assc(**kwargs)
    ns2slim_assoc = get_slim_ns2assoc(ns2assoc, dag, slim_dag)

    gene2go = objanno.get_id2gos_nss(**kwargs)
    go2genes = objanno.get_id2gos_nss(go2geneids=True, **kwargs)
    return dag, objanno, ns2assoc, gene2go, go2genes, slim_dag, ns2slim_assoc


# @st.cache_resource
def ontology_enrichment(
    query_genes: list[str], bg_genes: list[str], **kwargs
) -> tuple[GOEnrichmentStudyNS, list[GOEnrichmentRecord]]:
    """Perform ontology enrichment analysis."""
    oea_obj = GOEnrichmentStudyNS(bg_genes, **kwargs)

    # Run enrichment analysis
    oea_results = oea_obj.run_study(query_genes, **kwargs)
    oea_results_sig = [
        r
        for r in oea_results
        if (r.p_fdr_bh < kwargs["alpha"]) and (r.enrichment == "e")
    ]

    return oea_obj, oea_results_sig


# @st.cache_data
def format_ontology_enrichment_results(
    label: str, oea_results_sig: list[GOEnrichmentRecord], **kwargs
) -> pd.DataFrame:
    """Format ontology enrichment results into a pandas DataFrame."""
    # itemid2name to covert itemid to name
    # transform the result to a dataframe
    try:
        oea_results_sig_prt = pd.DataFrame(get_goea_nts_prt(oea_results_sig, **kwargs))
    except Exception as e:
        oea_results_sig_prt = create_enrichment_dataframe(oea_results_sig, **kwargs)

    if oea_results_sig_prt.empty:
        return pd.DataFrame()
    else:
        # calculate the gene_ratio and term_coverage
        oea_results_sig_prt["gene_ratio"] = round(
            oea_results_sig_prt["study_count"] / oea_results_sig_prt["study_n"], 2
        )
        oea_results_sig_prt["term_coverage"] = round(
            oea_results_sig_prt["study_count"] / oea_results_sig_prt["pop_count"], 2
        )

        # sort gene items
        oea_results_sig_prt["study_items"] = oea_results_sig_prt["study_items"].apply(
            lambda x: ",".join(sorted(x.split(", ")))
        )
        oea_results_sig_prt["pop_items"] = oea_results_sig_prt["pop_items"].apply(
            lambda x: ",".join(sorted(x.split(", ")))
        )

        # keep the columns we need and reorder them
        # prt_columns = ["GO", "NS", "enrichment", "name", "ratio_in_study", "ratio_in_pop", "p_uncorrected", "depth", "study_count", "p_fdr_bh", "study_items", "pop_items", "study_n", "pop_count", "pop_n", "item_id", "namespace", "level", "is_obsolete", "alt_ids", "defn"]
        kept_columns = [
            "NS",
            "GO",
            "name",
            "level",
            "depth",
            "p_uncorrected",
            "p_fdr_bh",
            "gene_ratio",
            "term_coverage",
            "study_count",
            "study_n",
            "pop_count",
            "pop_n",
            "ratio_in_study",
            "ratio_in_pop",
            "study_items",
            "pop_items",
        ]
        if "defn" in oea_results_sig_prt.columns:
            kept_columns.append("defn")
        oea_results_sig_prt = (
            oea_results_sig_prt[kept_columns]
            .copy()
            .rename(
                columns={
                    "NS": "namespace",
                    "GO": "term_id",
                    "name": "term",
                    "p_fdr_bh": "p_fdr",
                }
            )
        )
        return oea_results_sig_prt


@st.cache_data
def ontology_enrichment_pipeline(
    ontology_data: GeneOntologyData | FYPOData | MONDOData,
    query_genes: list[str],
    bg_genes: list[str],
    load_kwargs: dict = {},
    enrichment_kwargs: dict = {},
    format_kwargs: dict = {},
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Perform ontology enrichment analysis pipeline."""

    dag, objanno, ns2assoc, gene2go, go2genes, slim_dag, ns2slim_assoc = (
        load_ontology_data(ontology_data, **load_kwargs)
    )

    enrichment_kwargs["godag"] = dag
    enrichment_kwargs["ns2assoc"] = ns2assoc

    slim_enrichment_kwargs = enrichment_kwargs.copy()
    slim_enrichment_kwargs["godag"] = slim_dag
    slim_enrichment_kwargs["ns2assoc"] = ns2slim_assoc["all_ancestors"]
    slim_enrichment_kwargs["propagate_counts"] = False

    oea_obj, oea_results_sig = ontology_enrichment(
        query_genes, bg_genes, **enrichment_kwargs
    )
    oea_obj_slim, oea_results_sig_slim = ontology_enrichment(
        query_genes, bg_genes, **slim_enrichment_kwargs
    )

    oea_results_sig_prt = format_ontology_enrichment_results(
        "full", oea_results_sig, **format_kwargs
    )
    oea_results_sig_prt_slim = format_ontology_enrichment_results(
        "slim", oea_results_sig_slim, **format_kwargs
    )

    return oea_results_sig_prt, oea_results_sig_prt_slim


# @st.cache_data
def stringdb_api_functions(
    output_format: STRING_OUTPUT_FORMAT = "xml",
    method: STRING_METHOD = "get_string_ids",
    params: dict = {},
):
    """Perform STRING API functions."""
    request_url = "/".join([STRING_API_URL, output_format, method])
    for attempt in range(MAX_RETRIES):
        try:
            response = requests.post(request_url, data=params)
            response.raise_for_status()
            if "Something went wrong!" in response.text:
                raise RequestException(
                    "communication_error: Something went wrong! -- please contact STRING maintainers if the issue persists"
                )
            break
        except RequestException as e:
            if attempt < MAX_RETRIES - 1:
                time.sleep(RETRY_DELAY)
            else:
                raise ValueError(
                    f"Error performing STRING API functions: {e} after {MAX_RETRIES} retries"
                )

    match output_format:
        case "xml":
            return pd.read_xml(StringIO(response.text))
        case "tsv":
            return pd.read_csv(StringIO(response.text), sep="\t")
        case "tsv-no-header":
            return pd.read_csv(StringIO(response.text), sep="\t", header=None)
        case "json":
            return pd.read_json(StringIO(response.text))
        case _:
            raise ValueError(f"Invalid output format: {output_format}")


# @st.cache_data
def format_string_enrichment_results(
    enrichment_df: pd.DataFrame, query_genes: list[str], background_genes: list[str]
) -> pd.DataFrame:
    """Format STRING enrichment results into a pandas DataFrame."""

    renamed_column_names = {
        "term": "term_id",
        "category": "namespace",
        "description": "term",
        "p_value": "p_uncorrected",
        "fdr": "p_fdr",
        "gene_ratio": "gene_ratio",
        "term_coverage": "term_coverage",
        "number_of_genes": "study_count",
        "study_n": "study_n",
        "number_of_genes_in_background": "pop_count",
        "pop_n": "pop_n",
        "ratio_in_study": "ratio_in_study",
        "ratio_in_pop": "ratio_in_pop",
        "preferredNames": "study_items",
    }
    enrichment_df = enrichment_df.rename(columns=renamed_column_names)
    enrichment_df["study_n"] = len(query_genes)
    enrichment_df["pop_n"] = len(background_genes)
    enrichment_df["ratio_in_study"] = enrichment_df.apply(
        lambda row: f"{row['study_count']}/{row['study_n']}", axis=1
    )
    enrichment_df["ratio_in_pop"] = enrichment_df.apply(
        lambda row: f"{row['pop_count']}/{row['pop_n']}", axis=1
    )
    enrichment_df["gene_ratio"] = round(
        enrichment_df["study_count"] / enrichment_df["study_n"], 2
    )
    enrichment_df["term_coverage"] = round(
        enrichment_df["study_count"] / enrichment_df["pop_count"], 2
    )
    enrichment_df = enrichment_df[list(renamed_column_names.values())]

    namespace_description = {
        "Process": "Biological Process (Gene Ontology)",
        "Component": "Cellular Component (Gene Ontology)",
        "Function": "Molecular Function (Gene Ontology)",
        "PMID": "Reference Publications (PubMed)",
        "NetworkNeighborAL": "Local Network Cluster (STRING)",
        "KEGG": "KEGG Pathways",
        "RCTM": "Reactome Pathways",
        "COMPARTMENTS": "Subcellular Localization (COMPARTMENTS)",
        "Keyword": "Annotated Keywords (UniProt)",
        "InterPro": "Protein Domains and Features (InterPro)",
        "SMART": "Protein Domains and Features (SMART)",
    }
    enrichment_df["namespace"] = enrichment_df["namespace"].map(namespace_description)

    namespace_order = list(namespace_description.values())
    enrichment_df["namespace_order"] = enrichment_df["namespace"].map({
        v: i for i, v in enumerate(namespace_order)
    })
    enrichment_df.sort_values(by="namespace_order", inplace=True)
    enrichment_df.drop("namespace_order", axis=1, inplace=True)

    return enrichment_df


@st.cache_data
def stringdb_enrichment(query_genes, bg_genes):
    """Perform STRING enrichment analysis using the STRING API."""

    get_string_id_params = {
        "identifiers": STRING_SEPARATOR.join(bg_genes),
        "species": STRING_SPECIES_ID,
        "limit": 1,
        "echo_query": 1,
        "caller_identity": STRING_CALLER_IDENTITY,
    }

    get_string_id_response = stringdb_api_functions(
        output_format="xml", method="get_string_ids", params=get_string_id_params
    )
    get_string_id_response = get_string_id_response["stringId"].tolist()

    enrichment_params = {
        "identifiers": STRING_SEPARATOR.join(query_genes),
        "species": STRING_SPECIES_ID,
        "background_string_identifiers": STRING_SEPARATOR.join(get_string_id_response),
        "caller_identity": STRING_CALLER_IDENTITY,
    }
    # Parse results
    enrichment_df = stringdb_api_functions(
        output_format="xml", method="enrichment", params=enrichment_params
    )

    # format the results
    enrichment_df = format_string_enrichment_results(
        enrichment_df, query_genes, bg_genes
    )

    return enrichment_df


# @st.cache_data
def display_enrichment_results(enrichment_results: pd.DataFrame) -> alt.Chart:
    """Display enrichment results."""
    charts = []
    for ns, ns_results in enrichment_results.groupby("namespace", sort=False):
        chart = (
            alt.Chart(ns_results)
            .mark_circle()
            .encode(
                alt.X("gene_ratio:Q", axis=alt.Axis(grid=True), title="Gene ratio"),
                alt.Y(
                    "term:N",
                    axis=alt.Axis(
                        grid=True, labelLimit=500, title="Term", orient="right"
                    ),
                    sort=alt.EncodingSortField(field="gene_ratio", order="descending"),
                ),
                size=alt.Size("term_coverage:Q", title="Term coverage"),
                color=alt.Color("p_fdr:Q", title="FDR").scale(
                    scheme="yelloworangered", reverse=True
                ),
                tooltip=ns_results.columns.tolist(),
            )
        ).properties(
            title=f"Enrichment results for {ns}",
        )
        charts.append(chart)

    return alt.vconcat(*charts)
