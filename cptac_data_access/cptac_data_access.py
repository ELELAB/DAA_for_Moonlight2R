import cptac
import pandas as pd
import argparse
import sys
import os
from sklearn.decomposition import PCA
import seaborn as sns
import matplotlib.pyplot as plt 




parser = argparse.ArgumentParser()

parser.add_argument("-c", "--cancertype",
                    type=str,
                    required=True,
                    help="Cancer type short name")
parser.add_argument("-s", "--source",
                    type=str,
                    required=True,
                    help="Pipeline data source short name")
parser.add_argument("-p", "--plot",
                    action="store_true",
                    help="if set, data plots will be created")


def plot_distribution_boxplot(abundance_matrix, cancer_name, source, output_dir):
    """
    This function creates a multi-boxplot chart of protein abundance 
    for each sample in the dataset and save it as .png figure 

    Parameters:
    abundance_matrix (DataFrame) : abundance matrix (rows=proteins; cols=samples)
    cancer_name (str) : cancer type short name
    
    """

    plt.figure(figsize=(16, 6))
    sns.boxplot(data=abundance_matrix, orient='v', width=0.6, fliersize=1)
    plt.xticks(rotation=90, ha='right')
    plt.tick_params(axis='x', labelsize=5)
    plt.xlabel("Sample")
    plt.ylabel("Log2 Protein Abundance")
    plt.title("Boxplots of Protein Abundance per Sample")
    plt.tight_layout()

    output_path = os.path.join(output_dir, f'protein_abundance_boxplots_{cancer_name}_{source}.png')

    plt.savefig(output_path, dpi=300)
    plt.figure(figsize=(12, 6))

def plot_distribution_density(abundance_matrix, cancer_name, source, output_dir):
    """
    This function creates a distribtuon density plot of protein abundance 
    for each sample in the dataset and save it as .png figure 

    Parameters:
    abundance_matrix (DataFrame) : abundance matrix (rows=proteins; cols=samples)
    cancer_name (str) : cancer type short name
    
    """
    
    # remove unused indexes of the df 
    df = abundance_matrix.reset_index(drop=True)

    for sample in df.columns:
        sns.kdeplot(abundance_matrix[sample], label=sample, alpha=0.5)

    plt.title("Density Plot of Protein Expression per Sample")
    plt.xlabel("Log2 Protein Abundance")
    plt.ylabel("Density")
    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left',ncol=4, borderaxespad=0.,fontsize=5)
    plt.tight_layout()

    output_path = os.path.join(output_dir, f'protein_abundance_density_plots_{cancer_name}_{source}.png')

    plt.savefig(output_path, dpi=300)
    plt.figure(figsize=(12, 6))

def plot_pca(abundance_matrix, cancer_name, source, output_dir):
    """
    This function perform Principal Component Analysis (PCA) 
    on the protein abundance dataset and a PCA plot saved as .png 

    Parameters:
    abundance_matrix (DataFrame) : abundance matrix (rows=proteins; cols=samples)
    cancer_name (str) : cancer type short name
    
    """ 
    # remove unused indexes of the df
    df = abundance_matrix.reset_index(drop=True)

    # create a dataframe that includes for each sample ID the condition label (Tumor/Normal)
    sample_ids = df.columns.tolist()
    label_dict = {}
    for s in sample_ids:
        if ".N" in s :
            label_dict[s] = "N"
        else:
            label_dict[s] = "T"
    
    columns = ["Sample_ID", "Labels"]
    df_labels = pd.DataFrame(label_dict.items(), columns = columns)
    
    # Perform PCA
    pca = PCA(n_components=2)
    df_filtered = df.T
    pca_result = pca.fit_transform(df_filtered)
    
    # store PCA results into a dataframe
    pca_df = pd.DataFrame(pca_result, columns=['PC1', 'PC2'], index = df_filtered.index)
    # assign the labels to each sample by mapping them
    df_labels_indexed = df_labels.set_index('Sample_ID')
    pca_df['group'] = pca_df.index.map(label_dict)
    # plot PCA 
    plt.figure(figsize=(8,6))
    sns.scatterplot(data=pca_df, x='PC1', y='PC2', hue='group', palette='Set1') 
    plt.title('PCA of Protein Abudance')
    plt.xlabel(f'PC1 ({pca.explained_variance_ratio_[0]*100:.1f}%)')
    plt.ylabel(f'PC2 ({pca.explained_variance_ratio_[1]*100:.1f}%)')
    plt.grid(True)
    
    output_path = os.path.join(output_dir, f'pca_{cancer_name}_{source}.png')

    plt.savefig(output_path, dpi=300)

def main():

    args = parser.parse_args()
    try:
        all_cancer_types = cptac.get_cancer_options()
    except Exception as e:
        print(f'Error accessing cancer options: {e}')
        sys.exit(1)
    
    try:
        cancer_df = all_cancer_types.loc[args.cancertype]
        cancer_source_row = cancer_df.loc[args.source]
    except KeyError as e:
        print(f'Invalid cancer type or pipeline: {e}')
        sys.exit(1)

    selected_cancer_source_has_proteomics = cancer_source_row.apply(lambda x: 'proteomics' in x).any()
    if not selected_cancer_source_has_proteomics:
        print(f'The {args.source} pipeline is not designed to process quantitative proteomic data. Please choose another source pipeline')
        sys.exit(1)
        
    # download proteomic dataset for the selected cancer type and pipeline
    cancer_type_name = args.cancertype.capitalize()
    try:
        func_cancer_type = getattr(cptac, cancer_type_name)
        cptac_cancer_type = func_cancer_type()
        proteomic_cancer_type_cptac = cptac_cancer_type.get_proteomics(args.source)
    except AttributeError:
        print(f'cptac does not have a class for cancer type {cancer_type_name}')
        sys.exit(1)
    except Exception as e:
        print(f'Error downloading data for {args.cancertype}: {e}')
        sys.exit(1)

   
    if proteomic_cancer_type_cptac.empty:
        print("Downloaded proteomic dataset is empty")
        syst.exit(1)
     # transpose the abundance table to get proteins as cols and samples as proteins
    try:
        proteomic_dataset_transpose = proteomic_cancer_type_cptac.T
        # drop rows (proteins) with at least one NaN value across samples
        proteomic_dataset_transpose_clean = proteomic_dataset_transpose.dropna(axis =0, how = 'any')
    except Exception as e:
        print(f'Error processing proteomic dataset: {e}')
    
    
    output_dir = 'cptac_dataset_results'
    if not os.path.exists(output_dir):
        os.mkdir(output_dir)
    
    output_path = os.path.join(output_dir, f'abundance_matrix_{args.cancertype}_{args.source}.csv')
    # save abundance matrix
    proteomic_dataset_transpose_clean.to_csv(output_path)

    if args.plot:
        # plot boxplots
        plot_distribution_boxplot(proteomic_dataset_transpose_clean, args.cancertype, args.source, output_dir )
        # plot density plots
        plot_distribution_density(proteomic_dataset_transpose_clean, args.cancertype, args.source, output_dir)
        # plot PCA
        plot_pca(proteomic_dataset_transpose_clean, args.cancertype, args.source, output_dir)

if __name__ == "__main__":
    try:
        main()
    except Exception as e:
        print(f"ERROR: {e}", file=sys.stderr)
        sys.exit(1)