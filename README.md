# DAA_for_Moonlight2R
Differential abundance analysis for Moonlight2R

## Introduction
This repository contains:
- A folder `cptac_data_access/` which includes a script for the access to quantitative proteomic data of tumor and normal samples of different cancer types from CPTAC Data Portal and for data visualization.

## Requirments
To reproduce the results follow these steps:
- Clone the github repository
```bash 
git clone https://github.com/ELELAB/DAA_for_Moonlight2R.git
```
- create a new environment as follows:
```bash 
module load conda
conda create --prefix ./DAA_env python=3.11 r-base=4.3.3
```
- activate the environment
```bash 
module load conda
conda activate DAA_env/
```

- to run the scripts make sure to install these packages in the environment as follows: 
```bash
    conda install pandas
    conda install numpy
    conda install matplotlib
    conda install seaborn
    conda install scikit-learn
    pip install cptac
```

