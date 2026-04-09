# DAA_for_Moonlight2R
Differential abundance analysis (DAA) for Moonlight2R

## Introduction
This repository contains:
- A folder `cptac_data_access/` which includes a script for the access to quantitative proteomic data of tumor and normal samples of different cancer types from CPTAC Data Portal and for data visualization.
- A folder `DAA_limma/` which includes a script to perform (DAA) of protein using `limma` statistical method.
More details are provided within each subfolder.
- Two environment configuartion files, namely `environment-python.yml` and `environment-r.yml` that should be used for the setup of two virtual environments to run the Python and R scripts included in this repository.

## Requirments
To reproduce the results follow these steps:
- Clone the github repository
```bash 
git clone https://github.com/ELELAB/DAA_for_Moonlight2R.git
```
- Create two new environments for Python and R separately- using the environment configuration files already provided in the repository
```bash 
module load conda
conda env create --prefix ./DAA_py_env -f environment-python.yml
conda env create --prefix ./DAA_r_env -f environment-r.yml
```
- Activate one environment at time to run the python script for proteomic data access [cptac_data_access.py](/cptac_data_access/cptac_data_access.py) and the R script for DAA analysis [limma.R](/DAA_limma/limma.R). 
You can check if the two have environments have been created successfully as follows:
```bash 
module load conda
conda activate DAA_py_env/

conda deactivate
conda activate DAA_r_env/
```
- Navigate to the directory [cptac_data_access/](/cptac_data_access/) to perform proteomic data access
```bash 
cd cptac_data_access/
```
- Navigate to the directory [DAA_limma/](/DAA_limma/) to perform DAA
```bash 
cd DAA_limma/
```