"""
This script prepares the data for the gene page.
"""

# ================================= Imports =================================
import pandas as pd
import streamlit as st
from typing import List
from pathlib import Path
from goatools.obo_parser import GODag
from datetime import date

# ================================= Functions =================================

def get_gene_list_from_query(text: str) -> List[str]:
    """Parse gene list from text input."""
    if not text.strip():
        return []
    
    # Split by commas, newlines, or spaces
    genes = []
    for separator in [',', '\n', ' ']:
        if separator in text:
            genes = [g.strip() for g in text.split(separator) if g.strip()]
            break
    
    return genes

def validate_gene_ids(genes: List[str], gene_info: pd.DataFrame, gene_level_LFCs: pd.DataFrame) -> tuple[List[str], List[str], List[str], List[str], List[str]]:
    """Validate gene list against available data."""

    name2sysID = dict(zip(gene_info["gene_name"], gene_info["gene_systematic_id"]))
    covered_names = gene_level_LFCs["Name"].tolist()
    covered_sysIDs = gene_level_LFCs.index.tolist()

    valid_genes = []
    invalid_genes = []
    covered_genes = []
    covered_gene_sysIDs = []
    uncovered_genes = []

    # if the gene is a name, convert it to sysID
    for gene in genes:
        if gene in name2sysID.keys():
            valid_genes.append(gene)
            if gene in covered_names:
                covered_genes.append(gene)
                covered_gene_sysIDs.append(name2sysID[gene])
            else:
                uncovered_genes.append(gene)
        elif gene in name2sysID.values():
            valid_genes.append(gene)
            if gene in covered_sysIDs:
                covered_genes.append(gene)
                covered_gene_sysIDs.append(gene)
            else:
                uncovered_genes.append(gene)
        else:
            invalid_genes.append(gene)
    return valid_genes, invalid_genes, covered_genes, covered_gene_sysIDs, uncovered_genes

def sidebar_gene_input(gene_info_with_essentiality: pd.DataFrame, gene_level_LFCs: pd.DataFrame) -> tuple[List[str] | None, bool]:
    """Set the sidebar for the plot page."""
    input_form = st.sidebar.form("gene_input", clear_on_submit=False, border=True)
    input_form.subheader("Enter query genes:")
    gene_input = input_form.text_area("(comma or newline separated)", value="SPAC1002.09c\nSPAC3G9.12", height=300)
    submit_button = input_form.form_submit_button("Submit")
    if submit_button:
        gene_ids = get_gene_list_from_query(gene_input)
        valid_genes, invalid_genes, covered_genes, covered_gene_sysIDs, uncovered_genes = validate_gene_ids(gene_ids, gene_info_with_essentiality, gene_level_LFCs)
        input_form.badge(f"{len(gene_ids)} genes submitted", icon=":material/arrow_right_alt:", color="gray")
        input_form.badge(f"{len(valid_genes)} valid genes", icon=":material/check:", color="green")
        input_form.badge(f"{len(invalid_genes)} invalid genes", icon=":material/close:", color="red")
        input_form.badge(f"{len(covered_genes)} covered genes", icon=":material/check_circle:", color="blue")
        input_form.badge(f"{len(uncovered_genes)} uncovered genes", icon=":material/error:", color="orange")

        if invalid_genes:
            with st.sidebar.expander("Invalid genes", expanded=False):
                st.text("\n".join(invalid_genes))
    
        if uncovered_genes:
            with st.sidebar.expander("Uncovered genes", expanded=False):
                st.text("\n".join(uncovered_genes))
    else:
        covered_gene_sysIDs = None

    return covered_gene_sysIDs, submit_button

def assign_term_name(term_ID, term_dag):
    if term_ID in term_dag:
        return term_dag[term_ID].name
    else:
        return "No record for {}".format(term_ID)

def format_phaf_file(fypo_obo_file: Path, phaf_file: Path) -> Path:
    """Format the phaf file to the go style gaf file."""
    phaf_dag = GODag(str(fypo_obo_file))
    phaf = pd.read_csv(
        phaf_file, sep="\t"
    ).query(
        "(`Allele type` == 'deletion' or `Allele type` == 'disruption') and Condition.str.contains('FYECO:0000005')"
    )
    phaf["DB"] = "PomBase"
    phaf["DB_Object_ID"] = phaf["Gene systematic ID"]
    phaf["DB_Object_Symbol"] = phaf["Gene symbol"]
    phaf["Qualifier"] = ""
    phaf["GO_ID"] = phaf["FYPO ID"]
    phaf["DB:Reference"] = phaf["Reference"]
    phaf["Evidence"] = phaf["Evidence"]
    phaf["With"] = ""
    phaf["Aspect"] = "FYPO"
    phaf["DB_Object_Name"] = phaf["FYPO ID"].apply(assign_term_name, term_dag=phaf_dag)
    phaf["Synonym"] = ""
    phaf["DB_Object_Type"] = "protein"
    phaf["Taxon"] = "taxon:4896"
    phaf["Date"] = phaf["Date"].str.replace("-", "")
    phaf["Assigned_By"] = phaf["#Database name"]
    phaf["Annotation_Extension"] = phaf["Extension"]
    phaf["Gene_Product_Form_ID"] = ""
    reformat_phaf = phaf[["DB", "DB_Object_ID", "DB_Object_Symbol", "Qualifier", "GO_ID", "DB:Reference", "Evidence", "With", "Aspect",
             "DB_Object_Name", "Synonym", "DB_Object_Type", "Taxon", "Date", "Assigned_By", "Annotation_Extension", "Gene_Product_Form_ID"]].copy()

    with open(phaf_file.parent / "phaf_go_style.tsv", "w") as f:
        f.write(f"!gaf-version: 2.2\n!generated-by: Yusheng Yang\n!date-generated: {date.today().strftime('%Y-%m-%d')}\n!URL: https://www.pombase.org/monthly_releases/2025/pombase-2025-09-01/phenotypes_and_genotypes/pombase_phenotype_annotation.phaf.tsv\n!contact: yangyusheng@nibs.ac.cn\n")
    
    reformat_phaf.to_csv(phaf_file.parent / "phaf_go_style.tsv", sep="\t", index=False, header=False, mode="a")

    return phaf_file.parent / "phaf_go_style.tsv"

def format_mondo_gaf_file(mondo_obo_file: Path, mondo_gaf_file: Path) -> Path:
    """Format the mondo gaf file to the go style gaf file."""
    mondo_dag = GODag(str(mondo_obo_file))
    mondo = pd.read_csv(mondo_gaf_file, sep="\t")
    mondo["DB"] = "Pombase"
    mondo["DB_Object_ID"] = mondo["#gene_systematic_id"]
    mondo["DB_Object_Symbol"] = mondo["gene_name"]
    mondo["Qualifier"] = ""
    mondo["GO_ID"] = mondo["mondo_id"]
    mondo["DB:Reference"] = mondo["reference"]
    mondo["Evidence"] = ""
    mondo["With"] = ""
    mondo["Aspect"] = "MONDO"
    mondo["DB_Object_Name"] = mondo["mondo_id"].apply(assign_term_name, term_dag=mondo_dag)
    mondo["Synonym"] = ""
    mondo["DB_Object_Type"] = "protein"
    mondo["Taxon"] = "taxon:4896"
    mondo["Date"] = mondo["date"].fillna("2025-09-01").str.replace("-", "")
    mondo["Assigned_By"] = "PomBase"
    mondo["Annotation_Extension"] = ""
    mondo["Gene_Product_Form_ID"] = ""
    reformat_mondo = mondo[["DB", "DB_Object_ID", "DB_Object_Symbol", "Qualifier", "GO_ID", "DB:Reference", "Evidence", "With", "Aspect",
             "DB_Object_Name", "Synonym", "DB_Object_Type", "Taxon", "Date", "Assigned_By", "Annotation_Extension", "Gene_Product_Form_ID"]].copy()

    with open(mondo_gaf_file.parent / "mondo_go_style.tsv", "w") as f:
        f.write(f"!gaf-version: 2.2\n!generated-by: Yusheng Yang\n!date-generated: {date.today().strftime('%Y-%m-%d')}\n!URL: https://www.pombase.org/monthly_releases/2025/pombase-2025-08-01/ontologies_and_associations/human_disease_association.tsv\n!contact: yangyusheng@nibs.ac.cn\n")
    
    reformat_mondo.to_csv(mondo_gaf_file.parent / "mondo_go_style.tsv", sep="\t", index=False, header=False, mode="a")

    return mondo_gaf_file.parent / "mondo_go_style.tsv"