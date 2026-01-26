# Differential Abundance Analysis (DAA)

## Introduction
This folder contains script to perform Differential Abundance Analysis (DAA) of quantitative proteomic data across two sample group conditions. The script requires proteomic datasets already preprocessed (normalised and log2-transformed). 
DAA is used to identify and measure differences in protein levels between two or more groups such as Tumor vs Normal tissues or samples. The process usually involves the application of statistical methods to identify significant differences with associated p-values and fold changes. 

## DAA with limma
`limma` [(DOI: 10.18129/B9.bioc.limma)](https://bioconductor.org/packages/release/bioc/html/limma.html) is an R/Bioconductor package widely used for the analysis of -omics datasets and include models that allows differential abundance to be assessed. Over the past decades, `limma` has been used for gene discovery through differential expression analysis of microarray and high-throughput PCR data. 
Differential Abundance analysis for proteomics data is crucial for accurate detection of phenotype-specific proteins, which can be useful in biomedical applications such as biomarker and drug target discovery. `limma`is widely used for proteomics and a recent benchmarking study on the optimization of DAA tools for proteomics data suggests to use limma for DAA downstream analysis of different MS-detected proteomic datasets (performing well in any quantification settings). 

`limma` implements empirical Bayes moderated t-test and works well on relatively small sample sizes. Different experiment configurations can be designed based on the biological assumptions and questions of the investigation context. The most common setup is having samples from two different groups (such as case and control or tumor and normal). Samples are usually assumed to be independet but can be paired. The most common setup is to proceed with unpaired sample design.

### Limma steps
1. Load the protein abundance matrix as input, containing proteins as rows and samples as columns.
2. If the input abundance matrix comes along with a column containing the Ensembl Protein IDs for each protein entry, a filtering step is performed in order to remove duplicate gene names and keep only the canonical isoform associated to a specific Ensembl ID (according to Uniprot annotation) for each gene/protein entry. 
3. Design matrix creation: it represents the design of the experiment.
4. Linear Model Fit: limma fits a linear model for each protein to estimate the mean abundance in tumor and normal groups. The method used is the Mean Model for categorical vairables (based on least squares approach). It takes as input a matrix-like data object contaiing log-ratios or log-expression values (i.e. the abundance matrix); a design matrix with rows correspnding to samples and columns to coefficients to be estimated (e.g. 1/0 for tumor/normal). The function fits multiple linear models (one for each protein in the dataset). The model estimates the group means for each protein.
5. Contrast Matrix: after linear model fit a contrast matrix is used to define the way of mean comparison (e.g. tumorMean-normalMean).
6. Emprical Bayes Statistics: given a linear model fit, a moderated t-statistics is computed using the empirical bayes moderation of the standard error towards a global value. This method estimates the variance for each protein; a common variance trend across all proteins; shrinks each protein variance estimate toward this common trend. This method improve stability and accuracy of the statistical inference, reducing the impact of noise and instability of individual estimates. This method is useful when the number of features (e.g. proteins) in a dataset is much larger than the number of samples. This steps ensure that the variance of each protein is adjusted to avoid false positives and false negatives. This mainly affect p-values estimates.
7. The output is a table containing for each protein in the dataset, a log2FC, average expression and p-values. 

## Reproduce results
Activate the conda environment 
```bash 
conda activate DAA_env/
conda config --add channels conda-forge
conda config --add channels bioconda
conda install bioconductor-limma
conda install -c conda-forge r-optparse
conda install -c conda-forge r-dplyr
conda install -c conda-forge r-tidyverse
```
Run the analysis 
```bash 
Rscript limma.R -m "abundance_matrix.csv" -s "sample_labels.csv" -c "cancer type short name"
```
## Input preparation
If you have downloaded the data from the CPTAC Data Portal through the pyhton API and using the script [cptac_data_access.py](../cptac_data_access/cptac_data_access.py), you should be able to run the Differential Abundance Analysis using the files produced by the script. 

To make the script running on your own data instead, we suggest to have an input abundance matrix where:
- The first column should be named `Name` and should contain the gene names. 
- The shape of the matrix should be: nrows = # proteins and ncols = # samples.
- If an additional column containing Ensembl Protein IDs for each gene name is available to you- the column should be named `Database_ID`.
We also suggest that the input `sample_labels.csv` file should contain information as follows:
- The first column should be named `Sample_ID` and should contain the sample IDs.
- The second column should be named `Label_Group` and should contain the condition/group label assigned to each sample (e.g. Tumor or Normal)

## Working example 
An example is provided here to run Differential Abundance Analysis (DAA) for Colon Adenocarcinoma (coad) proteomic dataset.
Note that to run the following command you need to specify the actual location of the input files. Here, we assume that the abundance matrix 
and the sample label files are both stored in the default folder- created when running the [cptac_data_access.py](../cptac_data_access/cptac_data_access.py)

```bash 
Rscript limma.R -m ../cptac_data_access/cptac_dataset_results/abundance_matrix_coad_umich.csv -s ../cptac_data_access/cptac_dataset_results/sample_labels_coad_umich.csv -c coad
```
