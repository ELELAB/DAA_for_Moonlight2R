# CPTAC data access and proteomic dataset visual exploration

## CPTAC data access
This folder contains the `cptac_data_access.py` script which perform programmatic data access of quantitative proteomic dataset based on user-defined cancer type and pipeline data source.
Before running the script- look at the `cptac_cancer_info_table.csv` to get insights on the namings of each specific cancer type. In particular, you will need to choose- in the "Cancer Type Abbreviation" column- the name of the cancer type of interest- and in the "Sources" column- one of the sources pipelines that processed the -omics data. Note pipelines are tailored to the -omics data type and mainly "umich" and "bcm" are designed for proteomic data processing
For details on each pipeline please read the supplemental information in https://doi.org/10.1016/j.ccell.2023.06.009.

The `cptac_data_access.py` script will take as input the cancer type name and the source name and it will retrieve the corresponding proteomic abundance table, if available.

### Reproduce results
To reproduce the results follow these steps:
1) Activate environment 
```bash 
module load conda
conda activate "custom_env_name"
```
2) Check the  `cptac_cancer_info_table.csv` and choose one of the listed Cancer Type, its corresponding short name and the source pipeline name.
3) Run the script `cptac_data_access.py` as follows:

```python 
python cptac_data_access.py -c "cancer type short name" -s "pipeline source short name" 
```

## Visualization support
The `cptac_data_access.py` script also perform Principal Component Analysis (PCA) on the input dataset and create three figures:
- Boxplot distribution of protein abundance values of each sample across the entire dataset (protein_abundance_boxplots.png)
- Density plot of the proteins abundance values of each samples across the entire dataset (protein_abundance_density_plots.png)
- PCA plot where each sample is clustered on the group of belonging (e.g. Tumor, Normal) 

### Reproduce results
To reproduce the results with the plotting feature, run this command:

```python 
python cptac_data_access.py -c "cancer type short name" -s "pipeline source short name" -p
```

### Working example
To download and explore the Colon Adenocarcinoma proteomics dataset processed by the University of Michigan (UMICH) pipeline and get quantitative proteomics data- you can run the following:
```python 
python cptac_data_access.py -c coad -s umich -p
```
To download and explore the Colon Adenocarcinoma proteomics dataset processed by the Baylor College of Medicine (BCM) pipeline and get pan-cancer harmonized quantitative proteomic data- you can run the following:
```python 
python cptac_data_access.py -c coad -s bcm -p
```