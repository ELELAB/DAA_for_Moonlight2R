import cptac
import pandas as pd
import argparse
import sys
import os
from sklearn.decomposition import PCA
import seaborn as sns
import matplotlib.pyplot as plt 
import time  
import requests

parser = argparse.ArgumentParser()

parser.add_argument("-c", "--cancertype",
                    type=str,
                    required=True,
                    help="Cancer type short name.")
parser.add_argument("-s", "--source",
                    type=str,
                    required=True,
                    help="Pipeline data source short name.")
parser.add_argument("-p", "--plot",
                    action="store_true",
                    help="if set, data plots will be created.")
parser.add_argument("-i", "--isoform-mode",
                    type=str,
                    choices=["conservative", "strict"],
                    default="conservative",
                    help="Canonical isoform filtering. strict: keep a protein only if a canonical "
                         "isoform ENSP is in the data. conservative: if no canonical ENSP is in the "
                         "data, keep the row with the fewest missing values. Default: conservative.")


UNIPROT_SEARCH_URL = "https://rest.uniprot.org/uniprotkb/search"
GENE_NAME_ALIASES = ["name", "gene name", "gene names", "gene symbol", "symbol", "gene"]
DATABASE_ID_ALIASES = ["database id", "id", "ensembl id"]
CASE_DESCRIPTIONS = { "no_reviewed_entry": "the gene has no reviewed UniProt entry",
                      "not_in_uniprot": "none of the gene's ENSPs is linked to any isoform in UniProt",
                      "alternative_isoform": "all of the gene's ENSPs are linked to alternative (non-canonical) isoforms",
                      "mixed": "the gene's ENSPs are linked to alternative isoforms or not linked in UniProt; none is canonical",}

def find_column(df, aliases):
    """
    Return the name of the column matching one of the aliases, or None.
 
    Parameters:
    df (DataFrame) : input table
    aliases (list) : accepted names, in order of preference 
 
    Output:
    column name (str) or None if not found
    """
    # for each alias in the default alias list
    for alias in aliases:
        # for each column of the dataframe
        for col in df.columns:
            # fill any underscore or dot with empty spaces to parse the alias list
            simplified_str = " ".join(str(col).replace("_", " ").replace(".", " ").split()).lower()
            # if match is found, return the column name
            if simplified_str == alias:
                return col
    return None

def get_uniprot_data(session, url, params=None, timeout=120, retries=3, wait=10):
    """
    Send one request to UniProt. If the server has a temporary problem
    (5xx error, timeout, lost connection), wait and try again.
 
    Parameters:
    session (requests.Session) : open connection to reuse
    url (str)      : address to request
    params (dict)  : query parameters (None if already in the url)
    timeout (int)  : seconds to wait for an answer
    retries (int)  : maximum number of attempts
    wait (int)     : seconds to wait between attempts
 
    Output:
    the server response
    """
    for attempt in range(1, retries + 1):
        try:
            res = session.get(url, params=params, timeout=timeout)
            if res.status_code >= 500 and attempt < retries:
                print(f"UniProt server error {res.status_code}, retrying in {wait} s ({attempt}/{retries})")
                time.sleep(wait)
                continue
            res.raise_for_status()
            return res
        except (requests.ConnectionError, requests.Timeout) as e:
            if attempt == retries:
                raise
            print(f"Connection problem ({e}), retrying in {wait} s ({attempt}/{retries})")
            time.sleep(wait)

def download_uniprot_human():
    """
    This function downloads all reviewed human UniProt entries (only gene
    names, Ensembl links and isoform information) from the latest UniProt
    release, keeping the information in memory only. 
     
    Output:
    entries (list of dicts) : UniProt entries (dicts)
    """

    print("Downloading reviewed human UniProt entries...")
    # Query Uniprot datasets in pages by defining a set of parameters for the search 
    page_size=500
    session = requests.Session()
    url = UNIPROT_SEARCH_URL
    # search only for reviewed and human proteins, their accession ID and gene name.
    # Include ensembl cross references and alternative products annotations
    params = { "query": "reviewed:true AND organism_id:9606",
               "format": "json",
               "fields": "accession,gene_primary,gene_synonym,xref_ensembl,cc_alternative_products",
               "size": page_size,}
    
   
    entries = []
    release = "unknown"
    # start the download loop. In the first iteration, url is the first page address
    while url:
        # ask uniprot for the current page and store the answer in res
        # the function sends one request and if something goes wrong, waits and tries again (max 3 times)
        res = get_uniprot_data(session, url, params)
        release = res.headers.get("x-uniprot-release", release)
        # update the entries list with the Uniprot query results
        entries.extend(res.json().get("results", []))
        # Update the url with the address of the next page
        url = res.links.get("next", {}).get("url")
        # the next-page address already contains all parameters
        params = None  

    print(f"Downloaded {len(entries)} UniProt entries (UniProt release {release})")
    
    return entries

def get_uniprot_entries_gene_names(entries):
    """
    This function collects the gene names (main names and alternative names)
    of all reviewed UniProt entries
 
    Parameters:
    entries (list) : UniProt entries
 
    Output:
    list of gene names
    """
    names = []
    for entry in entries:
        if "genes" not in entry:
            continue
        for gene in entry["genes"]:
            if "geneName" in gene:
                names.append(gene["geneName"]["value"].upper())
            # check for alternative gene names 
            if "synonyms" in gene:
                for synonym in gene["synonyms"]:
                    names.append(synonym["value"].upper())

    return names

def get_uniprot_canonical_isoform(uniprot_data):
    """
    This function retrieves the Uniprot canonical isoform accession IDs 

    Parameters:
    uniprot data: full uniprot entry

    Output:
    canonical isoform ID or None if not found 
    """
    for element in uniprot_data.get("comments", []):
        if element.get("commentType") != "ALTERNATIVE PRODUCTS":
            continue
        for isoform in element.get("isoforms", []):
            if isoform.get("isoformSequenceStatus") == "Displayed":
                ids = isoform.get("isoformIds", [])
                if ids:
                    return ids[0]
    return None

def get_alternative_isoforms(uniprot_data):
    """
    True if the UniProt entry lists at least one isoform other than the canonical
    (any isoform whose status is not 'Displayed').
    """
    for element in uniprot_data.get("comments", []):
        if element.get("commentType") != "ALTERNATIVE PRODUCTS":
            continue
        for isoform in element.get("isoforms", []):
            if isoform.get("isoformSequenceStatus") != "Displayed":
                return True
    return False

def get_ensp_id(xref):  
    """
    Get the Ensembl protein ID (ENSP) from a UniProt cross-reference.
    """
    for prop in xref.get("properties", []):
        if prop.get("key") == "ProteinId":
            return prop.get("value")
    return None 

def strip_ensp_version(ids):
    """
    Remove the Ensembl version number: 'ENSPXXXX.Y' -> 'ENSPXXXX'
    """
    return ids.astype("string").str.replace(r"\.\d+$", "", regex=True)

def build_ensp_lookup_from_uniprot(entries, ignore_version=True):
    """
    This function builds a lookup table reporting for each Ensembl protein ID (ENSP) reported in UniProt, 
    whether is linked to the canonical isoform of the protein.
    If ENSP has no explicit assignment to an isoform and one single canonical isoform is reported in Uniprot -> reported as linked to the canonical isoform
    If ENSP has no explicit assignment to an isoform and multiple isoforms are reported in Uniprot -> reported as not linked to the canonical isoform 
    
    It returns a lookup table for each ENSP its True/False canonical isoform association.
 
    Parameters:
    entries (list)        : UniProt entries
    ignore_version (bool) : True  -> IDs without version ('ENSPXXXX')
                            False -> IDs with version   ('ENSPXXXX.Y')
 
    Output:
    Series: index = ENSP, value = True (canonical) / False
    """
    # iterate over the entries retrieved from uniprot 
    rows = []
    for entry in entries:
        # get the canonical isoform uniprot ID (if single isoform, the accession number is used directly)
        canonical_isoform = get_uniprot_canonical_isoform(entry) or entry.get("primaryAccession")
        # check if the protein lists multiple isoforms (alternatives)
        has_alternative_isoform = get_alternative_isoforms(entry)
        
        # get the cross-ref data (ensembl IDs) and its associated isoform uniprot ID
        for cross_ref in entry.get("uniProtKBCrossReferences", []):
            if cross_ref.get("database") != "Ensembl":
                continue
            ensp = get_ensp_id(cross_ref)
            if ensp is None:
                continue
            # reads which Uniprot isoform the ensp belongs to
            isoform = cross_ref.get("isoformId")
            # if not isoform id is linked to the ensp:
            # Case 1: protein has only one isoform reported (canonical). This might be the reason of missing isoform id -> ensp assigned to canonical 
            # Case 2: protein has multiple alternative isoforms -> ensp not assigned to canonical
            if isoform is None:
                # if protein has alternative isoforms- the current ensp may be not linked to canonical
                if has_alternative_isoform:
                    continue
                # if protein has no alternative isoforms- the reason why ensp is missing isoform ID may be linked to the canonical isoform
                isoform = canonical_isoform

            # add the selected ensp linked to the canonical isoform ID
            rows.append({"ensp": ensp, "is_canonical": isoform == canonical_isoform})
    
    # build the look up table 
    lookup_table = pd.DataFrame(rows, columns=["ensp", "is_canonical"])
    # if no exact match of the version is required, the ensp version is removed
    if ignore_version:
        lookup_table["ensp"] = strip_ensp_version(lookup_table["ensp"])
    # if an ensp is linked to more than one entry, it counts as canonical if it is canonical in any
    return lookup_table.groupby("ensp")["is_canonical"].any()
    

def get_sample_group_labels_df(abundance_matrix):
    """
    This function creates a dataframe containing for each sample ID the group label it belongs to.
    This is done following the nomenclature used from cptac.

    Parameters:
    abundance_matrix (DataFrame) : abundance matrix (rows=proteins; cols=samples)

    Output:
    out_df (DataFrame): sample ids mapped with group label
    """
    df_cols = abundance_matrix.columns
    
    sample_groups_list = []
    for id in df_cols:
        sample_label_pair = ()
        if id.endswith(".N"):
            sample_label_pair = (id, "N")
        else:
            sample_label_pair = (id,"T")
        
        sample_groups_list.append(sample_label_pair)
    
    out_df = pd.DataFrame(sample_groups_list, columns=["Sample_ID", "Label_Group"])
    return out_df


def choose_best_ensp(rows):
    """
    Among the rows of one gene (multiple ENSP), return the label of the row with the most
    non-missing values (ties are solved keeping first ENSP in alphabetical order).
    """
    ranked = rows.sort_values(["_n_values", "Database_ID"], ascending=[False, True])
    return ranked.index[0]

def choose_ensp_for_protein(ensps, is_reviewed, mode):
    """
    Apply the isoform rules to the rows (ensp) of each gene.

    Parameters:
    ensps (DataFrame)   : the rows of one gene, with the columns _status
                         ("canonical", "alternative" or "not_found") and _n_values
    is_reviewed (bool) : the gene has a reviewed UniProt entry
    mode (str)         : "conservative" or "strict"

    Output:
    best_ensp : label of the row to keep, or None if the gene is removed
    case (str): the gene category of belonging based on its own ENSPs annotations in uniprot 
    """
    # rule 1: genes without a reviewed entry are always removed
    if not is_reviewed:
        return None, "no_reviewed_entry"

    # rule 2: if the gene has canonical rows, keep the best of them (both modes) based on the number of NA values
    canonical_ensps = ensps[ensps["_status"] == "canonical"]
    if len(canonical_ensps) > 0:
        best_ensp = choose_best_ensp(canonical_ensps)
        case = "canonical"
        return best_ensp, case

    # no canonical row: if the gene has no ensp matching a canonical isoform
    # check if the gene has at least one ensp associated to alternative isoforms (conservative mode will accept it)
    has_alternative = (ensps["_status"] == "alternative").any()
    # check if the gene has at least one ensp not found in uniprot (not in the lookup table) (conservative mode will accept it)
    has_not_found = (ensps["_status"] == "not_found").any()
    # 1. If it has alternative isoforms and not found ensps, it will get the best among them (least NA) (only conservative mode)
    if has_alternative and has_not_found:
        case = "mixed"
    # 2. If it has only alternative ensps,  it will get the best among them (least NA) (only conservative mode)
    elif has_alternative:
        case = "alternative_isoform"
    # 3. If it has only not found ensps,  it will get the best among them (least NA) (only conservative mode)
    else:
        case = "not_in_uniprot"
    # rule 3: strict mode removes the gene. Conservative mode keeps its best row (the least missing values)
    if mode == "strict":
        return None, case
    best_ensp = choose_best_ensp(ensps)
    return best_ensp, case

def map_uniprot_isoforms(abundance_matrix, mode="conservative"):
    """
    This function filters the abundance matrix to keep one Ensembl Protein ID per gene name,
    preferring the UniProt canonical isoform. 

    The filtering can be performed in two modes:
    1) Conservative: The ENSPXXXX.Y is used by ingoring the exact versioning (Y)
        - If one or more of the gene's ENSPs match the canonical isoform in UniProt,
          the canonical ENSP with the fewest missing values is kept and all other
          ENSPs of the gene are dropped.
        - If none of the gene's ENSPs matches the canonical isoform (they match
          alternative isoforms, or are not found in UniProt), the gene is kept with
          its ENSP with the fewest missing values. This ENSP is not confirmed as
          the canonical isoform.   
    2) Strict: The ENSPXXXX.Y is used by considering the exact versioning (Y)
        - If one or more of the gene's ENSPs match the canonical isoform in UniProt,
          the canonical ENSP with the fewest missing values is kept and all other
          ENSPs of the gene are dropped.
    
    In both modes, only reviewed entries in Uniprot are kept.
    If ties on missing values, the first ENSP in alphabetical order is kept.
 
    Parameters:
    abundance_matrix (DataFrame) : abundance matrix (rows=proteins; cols=samples)
    mode (str)       : "conservative" or "strict"
 
    Output:
    out_df (DataFrame) : filtered abundance matrix (rows=proteins; cols=samples)
    """
    if mode not in ("conservative", "strict"):
        raise ValueError(f"Unknown isoform mode '{mode}'. Use 'conservative' or 'strict'.")
    
    # if any index name is present in the input df, convert indexes into normal columns
    # if no names in the index, any numerical index will be dropped 
    index_names = [n for n in abundance_matrix.index.names if n is not None]
    if index_names:
        df = abundance_matrix.reset_index()
    else:
        df = abundance_matrix.reset_index(drop=True)

    # parse the gene name and database ID columns in the input dataframe
    gene_name_col = find_column(df, GENE_NAME_ALIASES)
    database_id_col = find_column(df, DATABASE_ID_ALIASES)
 
    # if the gene name column is not found, the program should stop
    if gene_name_col is None:
        raise ValueError("No gene name column found")
    # if the database id is not found but the gene name is there, the program will not perform any isoform mapping 
    if database_id_col is None: # TODO: simulate this case
        print("No Database ID column found. The analysis will be run without UniProt isoform mapping")
        return abundance_matrix
    
    # assign standard names to the gene name and database id columns
    df = df.rename(columns={gene_name_col: "Name", database_id_col: "Database_ID"})

    # Download results from uniprot for all proteins and store the release data as well
    entries = download_uniprot_human()

    # Collect the gene names (including secondary gene names from uniprot)  
    reviewed_gene_names = get_uniprot_entries_gene_names(entries)
    
    # Build a lookup table to check the annotation at the isoform level of each ensp for each uniprot entry.
    # For each ensp in the input data, annotate its canonical isoform annotation status  
    isoform_annotation = None
    # CONSERVATIVE MODE 
    if mode == "conservative":
        # Build a lookup table to store for all the ensembl protein IDs reported in Uniprot the association to the canonical isoform  (True or False)
        lookup_ensp = build_ensp_lookup_from_uniprot(entries, True)
        # check if the ensp matches the uniprot ensp (without the version)
        database_id_no_version = strip_ensp_version(df["Database_ID"])
        # map the ensp of the input data with the uniprot ensp reported in the lookup table from uniprot entries
        isoform_annotation = database_id_no_version.map(lookup_ensp)
    # STRICT MODE
    elif mode == "strict":
        # Build a lookup table to store for all the ensembl protein IDs reported in Uniprot with the association to the canonical isoform  (True or False)
        lookup_ensp = build_ensp_lookup_from_uniprot(entries, False)
        # map the ensp of the input data with the uniprot ensp reported in the lookup table from uniprot entries (exact match)
        isoform_annotation = df["Database_ID"].astype("string").map(lookup_ensp)
    
    # Label each row with how its ENSP was found in UniProt
    # and count its non-missing values (the samples are the numeric columns)
    is_canonical = (isoform_annotation == True).fillna(False).astype(bool)
    is_alternative = (isoform_annotation == False).fillna(False).astype(bool)
    
    # counts non missing values per each gene_ensp
    df["_n_values"] = df.select_dtypes(include="number").notna().sum(axis=1)
    # place default as not found and update the columns for those entries associated with canonical or alternative isoform
    df["_status"] = "not_found"
    df.loc[is_canonical, "_status"] = "canonical"
    df.loc[is_alternative, "_status"] = "alternative"
    
    # Go through the genes one by one and apply the rules to each gene's own rows
    reviewed_gene_names = set(reviewed_gene_names)  
    kept_ensps = []
    removed_by_case = {}
    for gene_name, ensps in df.groupby("Name"):
        # the gene is reviewed if its name is a reviewed name,
        # or if any of its ENSPs was found in the reviewed UniProt entries
        found_in_uniprot = (ensps["_status"] != "not_found").any()
        is_reviewed = str(gene_name).upper() in reviewed_gene_names or found_in_uniprot
        # apply the filtering rules to keep only one ensp for the current gene
        kept_ensp, case = choose_ensp_for_protein(ensps, is_reviewed, mode)
        if kept_ensp is not None:
            kept_ensps.append(kept_ensp)
        else:
            if case not in removed_by_case:
                removed_by_case[case] = []
            removed_by_case[case].append(gene_name)

    
    # The filtered table contains the kept rows, in their original order, without the helper columns
    out = df.loc[sorted(kept_ensps)].drop(columns=["_status", "_n_values"])
    
    print(f"Canonical isoform mapping ({mode} mode): {len(df)} rows -> {len(out)} rows")
    for case, gene_names in removed_by_case.items():
        print(f"  {len(gene_names)} genes removed: ({CASE_DESCRIPTIONS[case]})")
    
    # Put the original column names back and restore the original layout
    out = out.reset_index(drop=True)
    out = out.rename(columns={"Name": gene_name_col, "Database_ID": database_id_col})
    if index_names:
        out = out.set_index(index_names)
 
    return out
 

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
        sys.exit(1)
     # transpose the abundance table to get proteins as cols and samples as proteins
    try:
        proteomic_dataset_transpose = proteomic_cancer_type_cptac.T

        print(f"Total number of proteins in the dataset: {len(proteomic_dataset_transpose)}")
        # canonical isoform mapping 
        matrix_canonical_isoform = map_uniprot_isoforms(proteomic_dataset_transpose, mode=args.isoform_mode) 
        # drop rows (proteins) with at least one NaN value across samples
        proteomic_dataset_transpose_clean = matrix_canonical_isoform.dropna(axis =0, how = 'any')
    except Exception as e:
        print(f'Error processing proteomic dataset: {e}')
        sys.exit(1)
    
    output_dir = 'cptac_dataset_results'
    if not os.path.exists(output_dir):
        os.mkdir(output_dir)
    
    output_path_abmat = os.path.join(output_dir, f'abundance_matrix_{args.cancertype}_{args.source}.csv')
    # save abundance matrix
    proteomic_dataset_transpose_clean.to_csv(output_path_abmat)

    # create and save sample groups labels dataframe
    sample_groups_df = get_sample_group_labels_df(proteomic_dataset_transpose_clean)
    output_path_samplab = os.path.join(output_dir, f'sample_labels_{args.cancertype}_{args.source}.csv')
    sample_groups_df.to_csv(output_path_samplab)
    

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