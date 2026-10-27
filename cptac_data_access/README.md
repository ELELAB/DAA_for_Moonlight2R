# CPTAC data access and proteomic dataset visual exploration

## CPTAC data access
This folder contains the `cptac_data_access.py` script which perform programmatic data access of quantitative proteomic dataset based on user-defined cancer type and pipeline data source.
Before running the script- look at the [cptac_cancer_info_table.csv](cptac_cancer_info_table.csv) to get insights on the namings of each specific cancer type. In particular, you will need to choose- in the "Cancer Type Abbreviation" column- the name of the cancer type of interest- and in the "Sources" column- one of the sources pipelines that processed the -omics data. Note pipelines are tailored to the -omics data type and mainly "umich" and "bcm" are designed for proteomic data processing.
For details on each pipeline please read the supplemental information in https://doi.org/10.1016/j.ccell.2023.06.009.

The `cptac_data_access.py` script will take as input the cancer type name and the source name and it will retrieve the corresponding proteomic abundance table, if available. 

## Canonical isoform mapping
For each protein in the input dataset, the script keeps a single Ensembl protein ID (ENSP ID), preferring the one linked to the canonical isoform as annotated in UniProt: the input ENSP IDs are matched against the Ensembl cross-references of reviewed human UniProt entries, according to the latest UniProt release.

**Modes**: 
The script ca be run in two different modes, depending on how strict the isoform filtering must be:
- conservative mode: keep ENSP IDs associated to alternative/unassigned isoform ID, if no canonical is avaialable. 
- strict mode: keep the ENSP IDs only if associated to the canonical isoform.

**Output**: 
A new filtered version of the cptac abundance matrix, reporting for each gene name and the selected ENSP ID, the corresponding UniProt AC.

*Note (1): if the input dataset does not contain the ENSP ID column, the isoform mapping will not be performed.*

**Isoform filtering design**
The script perform protein isoform mapping and filtering according to the following logic:
**1. Status of each ENSP ID.**
Each row of the dataset is matched to UniProt by its ENSP ID and its gene name (checking also for alternative gene names reported in the entries as a safety check when gene names assigned in the input data are not the updated ones).
According to this mapping, each ENSP ID gets a "status":
- **canonical**: ENSP ID linked to the canonical isoform of the entry. If Uniprot gives no isoform ID and the protein has a single isoform reported, it is assumed canonical.
- **alternative**: ENSP ID linked to alternative isoform.
- **unassigned**: ENSP ID is not linked to any isoform ID and the Uniprot entry reports other ENSP-isoform links.
- **not found** : ENSP ID is not reported in the UniProt entry.

*Note (2): If an input gene name and its associated ENSP ID appear in more than one UniProt entry (e.g. when the same protein entity is encoded by different genes- all sharing the same gene name in the input dataset, but each with its own   UniProt entry), the entry in which the ENSP is linked to canonical isoform is preferred, then alternative, then unassigned*

*Note (3): if an input gene name has multiple associated ENSP IDs each matching a different Uniprot entry (Uniprot ACs), all rows will be kept for that gene*

**2. Rows (ENSP IDs) always removed:**
- rows without a gene name in the input dataset.
- rows whose ENSP ID is in UniProt only under a different gene name (gene name mismatch).
- genes whose name is not in any reviewed UniProt entry.
- rows whose ENSP ID is not found in UniProt.

**3.ENSP ID selection and UniProt AC assignment**:
For each gene in the input dataset, the ENSP ID linked to it are assessed based on the assigned "status" (explained above): 
- If the gene has one or more canonical ENSP IDs, the one with the fewest missing value is kept. (This applies to both `strict` and `conservative` modes)¨
- If the gene has no canonical ENSP ID:
  1) `conservative` mode (default): the gene is kept, with its alternative or "unassigned" ENSP ID with the fewest missing values; if it has none, the gene is removed.
  2) `strict` mode: the gene is removed.

*Note*: ENSP versions (ENSPXXXX.Y, where Y = version) are ignored by default. If option `-e` is used, the full ENSP ID must match in UniProt to be considered valid.


### Reproduce results
To reproduce the results follow these steps:
1) Activate environment 
```bash 
module load conda
conda activate ../DAA_py_env/
```
2) Check the `cptac_cancer_info_table.csv` and choose one of the listed Cancer Type, its corresponding short name and the source pipeline name.

3) Run the script `cptac_data_access.py` follwing this command line template:

```bash 
python cptac_data_access.py -c "cancer type short name" -s "pipeline source short name" 
```
By default, the script will run in `conservative mode`. If you wish to run in `strict`mode run the command by adding the `-i` command-line option:
```bash 
python cptac_data_access.py -c "cancer type short name" -s "pipeline source short name" -i "strict"
```
To reproduce the results with exact ENSP ID matching, follow this command:
```bash 
python cptac_data_access.py -c "cancer type short name" -s "pipeline source short name" -e
```
## Visualization support
The `cptac_data_access.py` script also performs Principal Component Analysis (PCA) on the input dataset and create three figures:
- Boxplot distribution of protein abundance values of each sample across the entire dataset (protein_abundance_boxplots.png)
- Density plot of the proteins abundance values of each samples across the entire dataset (protein_abundance_density_plots.png)
- PCA plot where each sample is clustered on the group of belonging (e.g. Tumor, Normal) 

### Working example
To download and explore the Lung Squamous Cell Carcinoma  (LSCC) proteomics dataset processed by the University of Michigan (UMICH) pipeline and get quantitative proteomics data, you can run the following command:
```bash 
python cptac_data_access.py -c lscc -s umich -i conservative -p
```
The output is stored in a folder named `cptac_dataset_results` containing:
- downloaded cptac abundance matrix (`.csv`)
- sample labels mapping (`.csv`)
- PCA plot (`.png`)
- Abundance values distribution/density plots (`.png`)