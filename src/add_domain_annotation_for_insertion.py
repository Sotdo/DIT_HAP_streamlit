import argparse
from pathlib import Path

import numpy as np
import pandas as pd

def main(args):

    domains = pd.read_csv(args.domain_file, header=[0], sep="\t")

    insertion_annotation_list = []

    for file in args.insertion_dir.glob("*annotated_insertions.csv"):

        insertion_annotation_list.append(pd.read_csv(file, index_col=[0,1,2,3], header=0))
    
    insertion_annotations = pd.concat(insertion_annotation_list, axis=0)

    insertion_annotations.drop_duplicates(inplace=True)

    insertion_annotations[["domain_id", "domain_residues"]] = insertion_annotations.apply(
        lambda row: assign_protein_domain(row, domains), axis=1, result_type="expand"
    )

    insertion_annotation_with_domain = {}
    for file in args.ddr_dir.glob("*.csv"):
        sample = file.name.split(".")[0]
        DDRs = pd.read_csv(file, header=0)

        insertion_annotation_with_domain[sample] = insertion_annotations.reset_index().merge(DDRs.reset_index(), on=["Systematic ID", "domain_id", "domain_residues"], how="left").set_index(["#Chr", "Coordinate", "Strand","Target"])
        DDR_ratio = insertion_annotation_with_domain[sample] .groupby("Systematic ID")["DR"].transform(lambda x: x / x.max()).rename("DDR_ratio")

        insertion_annotation_with_domain[sample]["DDR_ratio"] = DDR_ratio.round(3)

        insertion_annotation_with_domain[sample].to_csv(args.output_folder / f"{sample}_with_domain_insertion_annotation.csv")

def assign_protein_domain(row, domain):
    Gene = row["Systematic ID"]
    if Gene in domain["Systematic ID"].values:
        residue = int(row["Residue_affected"])
        domains = domain[domain["Systematic ID"] == Gene].copy()
        for idx, row in domains.iterrows():
            for iDomain in row.loc["domain_residues"].split(","):
                domain_start = int(iDomain.split("-")[0])
                domain_end = int(iDomain.split("-")[1])
                domain_id = row.loc["domain_id"]
                if residue >= domain_start and residue <= domain_end:
                    return domain_id, row.loc["domain_residues"]
        return np.nan, np.nan
    else:
        return np.nan, np.nan

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Add domain annotation for insertion.")
    parser.add_argument(
        "-i",
        "--insertion-dir",
        dest="insertion_dir",
        required=True,
        type=Path,
        help="Dir of insertion",
    )
    parser.add_argument(
        "-d",
        "--domain-file",
        dest="domain_file",
        required=True,
        type=Path,
        help="File of domain",
    )
    parser.add_argument(
        "-ddr",
        "--DDR-dir",
        dest="ddr_dir",
        required=True,
        type=Path,
        help="Dir of DDR",
    )
    parser.add_argument(
        "-o",
        "--output",
        dest="output_folder",
        required=True,
        type=Path,
        help="Output Folder",
    )

    args = parser.parse_args()

    main(args)