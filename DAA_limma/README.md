# Differential Abundance Analysis (DAA)

## Introduction
This folder contains script to perform Differential Abundance Analysis (DAA) of quantitative proteomic data across two sample group conditions. The script requires proteomic datasets already preprocessed (normalised and log-transformed). 
DAA is used to identify and measure differences in protein leveles between two or more groups such as Tumor vs Normal tissue or samples. The process usually involves the application of statistical methods to identify significant differences with associated p-values and fold changes. 

## DAA with limma
`limma` is an R/Bioconductor package widely used for the analysis of -omics datasets and include models that allows differential abundance to be assessed. Over the past decades, `limma` has been used for gene discovery through differential expression analysis of microarray and high-throughput PCR data. 
Differential Abundance analysis for proteomics data is crucial for accurate detection of phenotype-specific proteins, which can be useful in biomedical applications such as biomarker and drug target discovery. `limma`is widely used for proteomics and a recent benchmarking study on the optimization of DAA tools for proteomics data suggests to use limma for DAA downstream analysis of different MS-detected proteomic datasets (performing well in any quantification settings). 

`limma` implements empirical Bayes moderated t-test and works well on relatively small sample sizes. Different experiment configurations can be designed based on the biological assumptions and questions of the investigation context. The most common setup is having samples from two different groups (such as case and control; tumor and normal). Samples are usually assumed to be independet but ca be paired. The most common setup is to proceed with unpaired sample design.

### Limma steps
1. Define a protein abundance matrix as input containing proteins as rows and samples as columns
2. Design matrix creation: it represents the design of the experiment.
3. Linear Model Fit: limma fits a linear model for each protein to estimate the mean abundance in tumor and normal groups. The method used is the Mean Model for categorical vairables (based on least squares approach). It takes as input a matrix-like data object contaiing log-ratios or log-expression values; a design matrix with rows correspnding to samples and columns to coefficients to be estimated (e.g. 1/0 for tumor/normal). The function fits multiple linear models (one for each protein in the dataset). The model estimate the group means for each protein.
4. Contrast Matrix: after linear model fit a contrast matrix is used to define the way of mean comparison (e.g. tumorMean-normalMean).
5. Emprical Bayes Statistics: given a linear model fit, a moderated t-statistics is computed using the empirical bayes moderation of the standard error towards a gloval value. This method estimate the variance for each protein; a common variance trend across all proteins; shrinks each protein variance estimate toward this common trend. This method improve stability and accuracy of statistical inference., reducing the impact of noise and instability of individual estimates. This method is useful when the number of features (proteins) in a dataset is much larger than the number of samples. This steps ensure that the variance of each protein- which indicates how much data iagress with itself- is adjusted to avoid false positive and false negative. This mainly affect p-values estimates.
6. The output is a table containing for each protein in the dataset, a log2FC, average expression and p-values. 
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

## Working results 