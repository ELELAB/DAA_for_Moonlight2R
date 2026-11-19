# python script
import argparse
import os
import sys

import numpy as np
import pandas as pd
from pathlib import Path
from pylimma import contrasts_fit, e_bayes, lm_fit, make_contrasts, top_table
from matplotlib.figure import Figure

RESULTS_DIR = "limma_results"

def load_data(data_path):
    # function to read a csv file
    return pd.read_csv(data_path)

def prepare_limma_input(abundance_matrix):
    '''
    This function reads the input abundance table and extracts the 
    necessary information required to run limma

    Parameters:
    abundance_matrix (dataframe): input proteomics abundance table

    Output:
    limma_input_matrix (dataframe): clean abundance matrix
    '''
    id_cols = ["ID", "Name", "Database_ID", "UniProt_ID"]
    limma_input_matrix = abundance_matrix.drop(
        columns=[c for c in id_cols if c in abundance_matrix.columns])
    limma_input_matrix.index = abundance_matrix["ID"]
    return limma_input_matrix

def run_limma(matrix, groups, contrast_labels):
    '''
    This function runs limma statistical test

    Parameters:
    matrix (dataframe): proteomics data 
    groups (dataframe): sample labels and group condition table
    contrast_labels (list): list containing the contrast strings
    
    Output: table with limma results

    '''
    # get the contrast labels to identify the group conditions and estimate the difference
    contrast_lab_1, contrast_lab_2 = contrast_labels[0], contrast_labels[1]
    labels = groups["Label_Group"].astype(str).to_numpy()

    if not (contrast_lab_1 in labels and contrast_lab_2 in labels):
        sys.exit("The specified labels for the contrast matrix design do not correspond to the "
                 "actual labels assigned to the samples. Please check labels definition and "
                 "specify them accordingly")

    # labels are assigned to the samples by row order
    if len(labels) != matrix.shape[1]:
        sys.exit(f"The matrix has {matrix.shape[1]} samples but the labels file has "
                 f"{len(labels)} rows")
    
    # check that the input sample labels exactly match the specified ones in the contrast label definition
    # and no other label is reported
    diff_labels = sorted(set(labels) - {contrast_lab_1, contrast_lab_2})
    if diff_labels:
        sys.exit(f"Samples with labels {diff_labels} are not in the contrast: limma.py cannot "
                 "handle them either. Please remove them from the input")

    # --- design matrix 
    # create two column names, each referring to a group condition
    design_cols = [f"groups{contrast_lab_1}", f"groups{contrast_lab_2}"]
    # create the design matrix, storing for each sample its group of belonging
    design = pd.DataFrame({
        design_cols[0]: (labels == contrast_lab_1).astype(float),
        design_cols[1]: (labels == contrast_lab_2).astype(float),
    })
    
    # fit the model to compute for each protein the average abundance in the two group conditions
    fit = lm_fit(matrix, design)
    
    # --- create the contrast estimation (e.g. average Tumor - average Normal)
    contrast_exp = f"{design_cols[0]} - {design_cols[1]}"
    contrasts_matrix = make_contrasts(contrasts=contrast_exp, levels=design)
    
    # compute the difference between the two averages (logFC)
    fit2 = contrasts_fit(fit, contrasts=contrasts_matrix)
    
    # For each protein- test whether its logFC is significantly different from 0 (t-statistic and p-value)
    fit3 = e_bayes(fit2, robust=True)

    # organise the results into a table
    limma_results = top_table(fit3, coef=0, adjust_method="BH", number=matrix.shape[0])
    limma_results = limma_results.rename(columns={
        "log_fc": "logFC", "ave_expr": "AveExpr", "p_value": "P.Value",
        "adj_p_value": "adj.P.Val", "b": "B"})
    limma_results.index.name = "ID"
    limma_results = limma_results.reset_index()

    limma_results["Direction"] = np.where(limma_results["logFC"] > 0, "up", "down")
    limma_results.loc[limma_results["logFC"].isna(), "Direction"] = np.nan
    return limma_results

def plot_limma_results(limma_table, threshold_pval, cancer_type):
    '''
    This function plots DAA results in the shape of a Volcano plot

    Parameters:
    limma_table (dataframe): output dataframe from limma analysis
    threshold_pval (float): threshold on the significance levels 
    cancer_type (str): string labels of the cancer type used
    '''
    if cancer_type != "None":
        filename = os.path.join(RESULTS_DIR, f"limma_volcano_plot_{cancer_type}.png")
    else:
        filename = os.path.join(RESULTS_DIR, f"limma_volcano_plot.png")

    table = limma_table.dropna(subset=["logFC", "adj.P.Val"])
    # as xlim(-3.1, 3.1) in R: points outside the range are not drawn
    table = table[(table["logFC"] >= -3.1) & (table["logFC"] <= 3.1)]
    sig = table["adj.P.Val"] < threshold_pval
    up = table["Direction"] == "up"
    categories = [
        (~sig & ~up, "black", "Downregulated insignificant"),
        (~sig & up, "black", "Upregulated insignificant"),
        (sig & ~up, "deepskyblue", "Downregulated significant"),
        (sig & up, "red", "Upregulated significant"),
    ]

    fig = Figure(figsize=(4.8, 4.8))
    ax = fig.subplots()
    for mask, colour, label in categories:
        sub = table[mask]
        ax.scatter(sub["logFC"], -np.log10(sub["adj.P.Val"]), s=2, c=colour, label=label)
    ax.set_xlim(-3.1, 3.1)
    ax.set_xlabel("log2(FC)", fontsize=15)
    ax.set_ylabel("-log10(p-value)", fontsize=15)
    ax.tick_params(labelsize=12)
    ax.spines[["top", "right"]].set_visible(False)
    ax.legend(loc="upper left", frameon=False, fontsize=7)
    fig.tight_layout()
    fig.savefig(filename, dpi=100)

# parse arguments
parser = argparse.ArgumentParser()
parser.add_argument("-i", "--input_matrix", 
                    type=str,
                    required=True, 
                    help="A csv file containing the protein abundance matrix")
parser.add_argument("-s", "--sample_labels",
                    type=str,
                    required=True,
                    help="A csv file containing for each sample, its group/condition label")
parser.add_argument("-d", "--design_contrast",
                    type=str,
                    default="T-N",
                    help="A string indicating how to design the contrast in the differential "
                         "abundance analysis (e.g. if T-N is specified- the fold change for each "
                         "protein will be the result of Tumor condition - Normal condition). "
                         "Please remember to use the same labels definition as in the input "
                         "provided in --sample_labels")
parser.add_argument("-c", "--cancer_type",
                    type=str,
                    default="None", 
                    help="A string indicating the cancer type")
parser.add_argument("-p", "--plot",
                    action="store_true",
                    default=False,
                    help="If the flag is specified, plots will be generated")
parser.add_argument("-t", "--threshold_pval",
                    type=float,
                    default=0.05,
                    help="If the flag is specified, a custom p value threshold will be applied "
                         "to the volcano plot")


def main():

    args = parser.parse_args()

    # parse contrast matrix labels
    contrast_labels = args.design_contrast.split("-")
    # load abundance matrix
    data_matrix = load_data(args.input_matrix)
    # load sample labels dataframe
    sample_labels_df = load_data(args.sample_labels)

    # identifier used by limma for each row:
    # the UniProt AC if the input has a UniProt_ID column, otherwise the gene name
    if "UniProt_ID" in data_matrix.columns:
        id_col = "UniProt_ID"
    elif "Name" in data_matrix.columns:
        id_col = "Name"
        print("No UniProt_ID column found: gene names are used as identifiers")
    else:
        sys.exit("The input matrix must contain a UniProt_ID or a Name column")
    
    data_matrix["ID"] = data_matrix[id_col]

    # every row must have an identifier, and only once
    missing_id = data_matrix["ID"].isna() | (data_matrix["ID"].astype(str) == "")
    if missing_id.any():
        sys.exit(f"{missing_id.sum()} rows have no {id_col} - please check the input matrix")

    if data_matrix["ID"].duplicated().any():
        duplicated_ids = data_matrix.loc[data_matrix["ID"].duplicated(), "ID"].unique()
        sys.exit(f"The input matrix contains duplicate {id_col} values: "
                f"{', '.join(map(str, duplicated_ids))}")

    # the identifier columns present in the input, added back to the limma results at the end
    annotation_cols = [c for c in ["Name", "Database_ID", "UniProt_ID"] if c in data_matrix.columns]
    annotation = data_matrix[["ID"] + annotation_cols]

    # run limma
    limma_input_matrix = prepare_limma_input(data_matrix)
    limma_results = run_limma(limma_input_matrix, sample_labels_df, contrast_labels)

    # remove Avg expression and B statistic columns from the output df (t is kept)
    limma_results = limma_results.drop(columns=["AveExpr", "B"])
    limma_results = limma_results.rename(columns={"P.Value": "PValue"})

    # add the identifier columns back, as first columns, and remove the helper column ID
    limma_results = limma_results.merge(annotation, on="ID", how="left")
    limma_results = limma_results[annotation_cols +
                              [c for c in limma_results.columns if c not in annotation_cols + ["ID"]]]


    cancer_type = args.cancer_type
    os.makedirs(RESULTS_DIR, exist_ok=True)

    # as write.csv(row.names = TRUE) in R: first column with row numbers 1..n
    limma_results.index = np.arange(1, len(limma_results) + 1)
    if cancer_type != "None":
        out_file = os.path.join(RESULTS_DIR, f"limma_table_{cancer_type}.csv")
    else:
        out_file = os.path.join(RESULTS_DIR, f"limma_table.csv")
    limma_results.to_csv(out_file, index=True, index_label="", na_rep="NA")

    if args.plot:
        plot_limma_results(limma_results, args.threshold_pval, cancer_type)

if __name__ == "__main__":
    main()