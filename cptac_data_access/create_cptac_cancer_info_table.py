# this script will create a .csv file which includes- for each cancer type- a list of pipeline data sources for proteomic data processing
import cptac
import pandas as pd


out_df = pd.DataFrame()
cancer_names = []
cancer_abbreviation = []
cancer_info = cptac.get_cancer_info()
for cancer_abbr in cancer_info.keys():
    if cancer_abbr == "pda":
        continue
    cancer_names.append(cancer_info[cancer_abbr])
    cancer_abbreviation.append(cancer_abbr)

out_df['Cancer Type Name'] = list(cancer_names)
out_df['Cancer Type Short Name'] = cancer_abbreviation
cancer_options = cptac.get_cancer_options()
cancer_options_abbrev = cancer_options.index.get_level_values(0).unique().drop('all_cancers')
cancer_data_sources = []
out_df['Sources'] = pd.Series(dtype='object')
i = 0
for cancer_type in cancer_options_abbrev:
    single_cancer_df = cancer_options.loc[cancer_type]
    cancer_data_sources = list(single_cancer_df.index)
    proteomics_cancer_data_sources = []
    for source in cancer_data_sources:
        if "proteomics" in single_cancer_df.loc[source][0]:
            proteomics_cancer_data_sources.append(source)
    out_df.at[i, 'Sources'] = proteomics_cancer_data_sources
    i+=1

out_df.to_csv('cptac_cancer_info_table.csv')