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
                         "data, keep the ensp linked to alternative or unassigned isoform. Default: conservative.")
parser.add_argument("-e", "--version-id",
                    action="store_true",
                    help="If set, the full version of the Database ID (ENSEMBL Protein ID) will be parsed for the isoform mapping")


UNIPROT_SEARCH_URL = "https://rest.uniprot.org/uniprotkb/search"
GENE_NAME_ALIASES = ["name", "gene name", "gene names", "gene symbol", "symbol", "gene"]
DATABASE_ID_ALIASES = ["database id", "id", "ensembl id"]
CASE_DESCRIPTIONS = { "no_reviewed_entry": "the gene has no reviewed UniProt entry",
                      "not_in_uniprot": "none of the gene's ENSPs is linked to any isoform in UniProt according to Ensembl cross-reference",
                      "alternative_isoform": "all of the gene's ENSPs are linked to alternative (non-canonical) or unassigned isoforms",}


def find_column(df, aliases):
    """
    This function returns the name of the column matching one of the aliases, or None.
 
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
    This function sends one request to UniProt. If the server has a temporary problem
    (e.g. timeout, lost connection), it waits and try again.
 
    Parameters:
    session        : open connection to reuse
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
               "fields": "accession,gene_primary,gene_synonym,gene_orf,xref_ensembl,cc_alternative_products",
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
    of all the input UniProt entries
 
    Parameters:
    entries (list) : UniProt entries
 
    Output:
    gene names (set)
    """
    names = set()
    for entry in entries:
        names |= get_entry_gene_names(entry)
    return names

def get_uniprot_canonical_isoform(uniprot_data):
    """
    This function retrieves the Uniprot canonical isoform accession IDs 

    Parameters:
    uniprot_data: full uniprot entry

    Output:
    canonical isoform ID (str) or None if not found  
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
    This function retrieves alternative isoform linked to the input UniProt entry.
    
    Parameter: Uniprot entry 

    Output:
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

def get_ensp_id(cross_ref):  
    """
    This function retrieves the Ensembl protein ID (ENSP) from a UniProt cross-reference.

    Output:
    ENSP ID (str)
    """
    for prop in cross_ref.get("properties", []):
        if prop.get("key") == "ProteinId":
            return prop.get("value")
    return None 

def strip_ensp_version(ids):
    """
    This function removes the Ensembl version number: 'ENSPXXXX.Y' -> 'ENSPXXXX'

    Parameter: ENSP ID

    Output:
    Stripped version of the ENSP ID (str)
    """
    return ids.astype("string").str.replace(r"\.\d+$", "", regex=True)

def get_entry_gene_names(entry):
    """
    This function retrieves all the gene names of one UniProt entry (current, alternative and ORF names)

    Parameters: UniProt entry

    Output:
    Gene names associated with the Uniprot entry (set)
    """
    names = set()
    for gene in entry.get("genes", []):
        if "geneName" in gene:
            names.add(gene["geneName"]["value"].upper())
        for synonym in gene.get("synonyms", []):
            names.add(synonym["value"].upper())
        for orf_name in gene.get("orfNames", []):
            names.add(orf_name["value"].upper())
    return names

def build_ensp_lookup_from_uniprot(entries, ignore_version=True):
    """
    This function builds a lookup table reporting for each Ensembl protein ID (ENSP ID)
    linked to a Uniprot entry, its Uniprot AC and its whether is associated to canonical or alternative isoform.
    
    For each entry-it reads:
    - accession: the Uniprot AC
    - is_canonical: whether the ENSP is canonical isoform of that entry -> True (linked to canonical); False (linked to alternative); None (listed in the entry, but not assigned to any isoform)
    - gene_names: the gene names of the entry (including alternative and ORF names)

    If the ENSP has no isoform ID in Uniprot:
    - the protein has a single isoform -> the ENSP assumed to be the canonical one (is_canonical = True)
    - the protein has several isoforms -> the isoform is unknown (is_canonical = None)
 
    If ignore_version is True, the version is removed from the ENSP IDs
    
    Note: If one single ENSP ID is associated to multiple entries- all the entries information will be reported for that ENSP ID
    
    Parameters:
    entries (list)        : UniProt entries
    ignore_version (bool) : True  -> IDs without version ('ENSPXXXX')
                            False -> IDs with version   ('ENSPXXXX.Y')

    Output:
    A dictionary: key = ENSP ID, value = list of the UniProt entries mentioning it,
        each as {"accession": str, "is_canonical": bool/None, "gene_names": set}
    """
    # iterate over the entries retrieved from uniprot 
    lookup = {}
    for entry in entries:
        # get the canonical isoform uniprot ID (if single isoform, the accession number is used directly)
        canonical_isoform = get_uniprot_canonical_isoform(entry) or entry.get("primaryAccession")
        entry_gene_names = get_entry_gene_names(entry)
        # check if the protein lists multiple isoforms (alternatives)
        has_alternative_isoform = get_alternative_isoforms(entry)
        
        # get the cross-ref data (ensembl IDs) and its associated isoform uniprot ID
        for cross_ref in entry.get("uniProtKBCrossReferences", []):
            if cross_ref.get("database") != "Ensembl":
                continue
            ensp = get_ensp_id(cross_ref)
            if ensp is None:
                continue
            # reads which Uniprot isoform the ensp belongs to (None if Uniprot does not specify it)
            isoform = cross_ref.get("isoformId")
            if isoform is not None:
                # ensp associated to an isoform: canonical or not
                is_canonical = (isoform == canonical_isoform)
            elif not has_alternative_isoform:
                # ensp not associated to an isoform ID and the enry has no alternative isoform -> single isoform: it is the canonical one
                is_canonical = True
            else:
                # ensp not associated to an isoform ID, and the entry has several isoforms -> isoform annotation in unknown
                is_canonical = None

            if ignore_version:
                ensp = ensp.split(".")[0]
            if ensp not in lookup:
                lookup[ensp] = []
            lookup[ensp].append({"accession": entry.get("primaryAccession"), "is_canonical": is_canonical, "gene_names": entry_gene_names })

    return lookup
    
def get_sample_group_labels_df(abundance_matrix):
    """
    This function creates a dataframe containing for each sample ID the group label it belongs to.
    This is done following the nomenclature used from cptac.

    Parameters:
    abundance_matrix (Dataframe) : abundance matrix (rows=proteins; cols=samples)

    Output:
    A Dataframe: sample ids mapped with group label
    """
    df_cols = abundance_matrix.columns
    
    sample_groups_list = []
    for id in df_cols:
        if id.endswith(".N"):
            sample_label_pair = (id, "N")
        else:
            sample_label_pair = (id,"T")
        
        sample_groups_list.append(sample_label_pair)
    
    out_df = pd.DataFrame(sample_groups_list, columns=["Sample_ID", "Label_Group"])
    return out_df

def choose_best_ensp(rows):
    """
    Among the given rows of one gene (multiple ENSP), return the label of the row with the most
    non-missing values (ties are solved by keeping first ENSP in alphabetical order).
    """
    ranked = rows.sort_values(["_n_values", "Database_ID"], ascending=[False, True])
    return ranked.index[0]

def choose_ensp_for_protein(ensps, is_reviewed, mode):
    """
    Apply the isoform rules to the rows (ensp) of a given gene:
    1. the gene has no reviewed UniProt entry -> removed (in both modes)
    2. the gene has one or more ENSP-canonical rows -> the best canonical row (with fewest missing values) of each Uniprot entry is kept (both modes)
    3. no canoncial rows -> the gene is removed (strict mode); the best (with fewest missing values) available row with a Uniprot AC is kept (alternative or unassigned)(conservative mode)

    Parameters:
    ensps (DataFrame)   : the rows of one gene, with the columns _status
                         ("canonical", "alternative", "not_found", "unassigned"), _n_values, _accession
    is_reviewed (bool) : the gene has a reviewed UniProt entry
    mode (str)         : "conservative" or "strict"

    Output:
    A collection (list) of labels (indexes) of the row to keep, or empty if the gene is removed
    case (str) : why: "no_reviewed_entry", "canonical", "multiple_proteins",
                  "alternative_isoform", "not_in_uniprot"
    """
    # rule 1: genes without a reviewed entry are always removed
    if not is_reviewed:
        return [], "no_reviewed_entry"

    # rule 2: if the gene has canonical rows, keep the best of them (both modes) based on the number of NA values
    canonical_ensps = ensps[ensps["_status"] == "canonical"]
    if len(canonical_ensps) > 0:
        kept = [choose_best_ensp(rows) for _, rows in canonical_ensps.groupby("_accession")]
        case = "canonical" if len(kept) == 1 else "multiple_proteins"
        return kept, case
    
    # rule 3:
    # no canonical row: if the gene has no ensp matching a canonical isoform
    # check if the gene has at least one ensp associated to alternative isoforms or unassigned isoforms (conservative mode accepts it)
    has_alternative = ensps["_status"].isin(["alternative", "unassigned"]).any()
    if has_alternative:
        case = "alternative_isoform" # at least one alternative or unassigned ENSP
    else:
        case = "not_in_uniprot" # no ENSP of the gene found in UniProt
    # strict mode removes the gene; conservative mode keeps its best row (the least missing values)
    if mode == "strict":
        return [], case

    # In conservative mode keep the best row among those with a UniProt AC (alternavie or unassigned)
    # rows whose ENSP is not in UniProt are never kept
    rows_with_ac = ensps[ensps["_status"].isin(["alternative", "unassigned"])]
    if len(rows_with_ac) == 0:
        return [], case
    return [choose_best_ensp(rows_with_ac)], case

def map_uniprot_isoforms(abundance_matrix, mode="conservative", use_version=False):
    """
    This function filters the abundance matrix to keep one Ensembl Protein ID per gene name,
    preferring the ENSP ID liked to the UniProt canonical isoform. 
    ENSP matching: by default, ENSP IDs are compared without their version (ENSPXXXX).
    If use_version=True, the version must match exactly (ENSPXXXX.Y). This applies in both modes.

    The filtering can be performed in two modes:
    1) Conservative:
        - If one or more of the gene's ENSPs match the canonical isoform in UniProt,
          the canonical ENSP with the fewest missing values is kept and all other
          ENSPs of the gene are dropped.
        - If none of the gene's ENSPs matches the canonical isoform, the alternative or unassigned ENSP
          with fewest missing values is kept; if there is none, the gene is removed. 
    2) Strict:
        - If one or more of the gene's ENSPs match the canonical isoform in UniProt,
          the canonical ENSP with the fewest missing values is kept and all other
          ENSPs of the gene are dropped.
        - If none of the ENSP are linked to a canonical isoform, the gene is removed.
    
    In both modes:
    - only reviewed entries in Uniprot are kept;
    - rows whose ENSP is in UniProt only under other gene names are removed, with a warning;
    - rows whose ENSP is not in UniProt are never kept, so every output row has a UniProt AC;
    - if a gene encodes several proteins with separate UniProt entries one canonical ENSP is kept for each; they share the gene name
      and are told apart by Database_ID and UniProt_AC;
    - if a row's ENSP and gene name match several UniProt entries (e.g. the same protein
      encoded by several genes), canonical entries are preferred and the first UniProt AC in alphabetical order is used, with a warning;
    - If for one gene in the dataset, multiple ENSP are associated to canonical isoforms of different UniProt entries (or alternative only) 
      the entry with fewest missing values is kept (If ties on missing values, the first ENSP in alphabetical order is kept)
 
    For each row (protein entry) the gene name, database ID and UniProt AC are kept in the final dataframe
    
    Parameters:
    abundance_matrix (DataFrame) : abundance matrix (rows=proteins; cols=samples)
    mode (str)       : "conservative" or "strict"
    use_version (bool): True -> ENSP IDs matched with their version; False -> without (default)
 
    Output:
    A DataFrame: filtered abundance matrix (rows=proteins; cols=samples)
    """
    if mode not in ("conservative", "strict"):
        raise ValueError(f"Unknown isoform mode '{mode}'. Use 'conservative' or 'strict'.")
    
    # if any index name is present in the input df, convert indexes into normal columns
    # if no names in the index, any numerical index will be dropped 
    index_names = [n for n in abundance_matrix.index.names if n is not None]
    if index_names:
        input_df = abundance_matrix.reset_index()
    else:
        input_df = abundance_matrix.reset_index(drop=True)

    # parse the gene name and database ID columns in the input dataframe
    gene_name_col = find_column(input_df, GENE_NAME_ALIASES)
    database_id_col = find_column(input_df, DATABASE_ID_ALIASES)
 
    # if the gene name column is not found, the program should stop
    if gene_name_col is None:
        raise ValueError("No gene name column found")
    # if the database id is not found but the gene name is there, the program will not perform any isoform mapping 
    if database_id_col is None: 
        print("No Database ID column found. The analysis will be run without UniProt isoform mapping")
        return abundance_matrix
    

    # assign standard names to the gene name and database id columns
    input_df = input_df.rename(columns={gene_name_col: "Name", database_id_col: "Database_ID"})
    # number of rows before any filtering (for the summary at the end)
    n_input_rows = len(input_df)
    n_input_genes = input_df["Name"].nunique()
    # rows without a gene name cannot be assigned to a protein: remove them
    has_name = input_df["Name"].notna()
    n_no_name = (~has_name).sum()
    input_df = input_df[has_name]
    
    # Download results from uniprot for all proteins and print the release version as well
    reviewed_entries = download_uniprot_human()
    # For all the entries- collect the gene names (including secondary gene names from uniprot)  
    reviewed_gene_names = get_uniprot_entries_gene_names(reviewed_entries)
    
    # Build a lookup table to store the annotation at the isoform level of each ensp for each uniprot entry.
    if use_version:
        # full ENSP ID, with version (e.g. ENSP00000269305.4)
        lookup_ensp = build_ensp_lookup_from_uniprot(reviewed_entries, ignore_version=False)
        input_ensp = input_df["Database_ID"].astype("string")
    else:
        # ENSP ID without version (e.g. ENSP00000269305)
        lookup_ensp = build_ensp_lookup_from_uniprot(reviewed_entries, ignore_version=True)
        input_ensp = strip_ensp_version(input_df["Database_ID"])     
    
    # add a column in the input dataframe to keep track of the count of non-missing value for each row
    input_df["_n_values"] = input_df.select_dtypes(include="number").notna().sum(axis=1)

    # for each row of the input dataset, find the Uniprot entry with the same ENSP ID and the same gene name in the lookup table
    # this step is done by adding two new temporary columns in the input dataframe: 
    # status: not found/canonical/alternative/unassigned
    # accession: Uniprot AC
    input_df["_status"] = "not_found"
    input_df["_accession"] = None
    # keep track of those cases in which the ENSP ID is found in Uniprot, but the gene name does not match
    gene_name_mismatch = []
    # list for warning when a row is kept only as a "backup" (alernative, unassigned)
    kept_non_canonical = []
    for label, ensp, gene_name in zip(input_df.index, input_ensp, input_df["Name"]): 
        # get the entries associated to the current ensp ID in the lookup table (there could be multiple UniProt entries associated)
        # if the entry is not found, it returns an empty list
        candidates = lookup_ensp.get(ensp,[]) 
        # among these entries, check which entries' gene names include the current ensp gene name
        gene_name_match = [c for c in candidates if str(gene_name).upper() in c["gene_names"]]
        # the ENSP is in Uniprot, but only under other gene names: these cases are dropped
        if len(candidates) > 0 and len(gene_name_match) == 0: # if candidates not empty (ENSP ID found in lookup) but gene_name_match empty (no gene name match)
            gene_name_mismatch.append(label)
            print(f"Warning: {gene_name} / {input_df.loc[label, 'Database_ID']} removed: gene name not found for this ENSP ID in Uniprot")
        # if the ENSP is in Uniprot under the current row gene name- two cases:
        # 1. one single match -> one protein -> its status is used (canonical/alternative)
        # 2. one ENSP ID matches to multiple Uniprot IDs and gene names (e.g. protein families) -> assign the uniprot AC of the entry that maps the ENSP ID as canonical.
        # 2a. if no canonical -> assingn the alternative isoform uniprot AC 
        # 2b. if no alternative -> assign the uniprot AC reporting the ENSP with no isoform annotation (unassigned)
        # 2c. If multiple matches exists - the first one is taken in alphabetical order
        elif len(gene_name_match) > 0:
            # among the matching entries, prefer canonical, then alternative, then unassigned
            canonical_matches = [c for c in gene_name_match if c["is_canonical"] is True]
            alternative_matches = [c for c in gene_name_match if c["is_canonical"] is False]
            if len(canonical_matches) > 0:
                # the ENSP is canonical in at least one matching entry
                selected_matches = canonical_matches
                input_df.loc[label, "_status"] = "canonical"
            elif len(alternative_matches) > 0:
                # the ENSP is an alternative isoform in the matching entries
                selected_matches = alternative_matches
                input_df.loc[label, "_status"] = "alternative"
            else:
                # the ENSP is listed in the entry, but not assigned to any isoform
                selected_matches = gene_name_match
                input_df.loc[label, "_status"] = "unassigned"

            # in the case multiple Uniprot IDs are associated to the same ENSP ID (e.g. protein families)
            # the first is kept in alphabetical order
            selected_accessions = sorted({c["accession"] for c in selected_matches})
            input_df.loc[label, "_accession"] = selected_accessions[0]
            # warn when several uniprot entries matched
            all_accessions = sorted({c["accession"] for c in gene_name_match})
            if len(all_accessions) > 1:
                print(f"Warning: {gene_name} / {input_df.loc[label, 'Database_ID']}: this ENSP and gene name "
                      f"match {len(all_accessions)} UniProt entries ({', '.join(all_accessions)}). "
                      f"{input_df.loc[label, '_status'].capitalize()} status and UniProt ID "
                      f"{selected_accessions[0]} are used.")
    
    # number of genes with at least one mismatched row 
    n_mismatch_genes = input_df.loc[gene_name_mismatch, "Name"].nunique()
    # remove those rows which had a gene name mismatch in the dataframe
    input_df = input_df.drop(index=gene_name_mismatch)
    
    # Go through the genes one by one and apply the rules to each gene's own rows 
    kept_ensps = []
    removed_genes_by_case = {}
    removed_rows_by_case = {} 
    n_extra_rows = 0 # number of rows removed because only one ensp (ie the canonical one) is kept
    n_extra_genes = 0   # kept genes that had other ENSPs dropped
    for gene_name, ensps in input_df.groupby("Name"):
        # the gene is reviewed if its name (main, alternative or ORF name) belongs to a reviewed UniProt entry
        is_reviewed = str(gene_name).upper() in reviewed_gene_names
        # apply the isoform filtering rules for the current gene
        # the function returns the list of row labels (indexes) and the case collection
        kept, case = choose_ensp_for_protein(ensps, is_reviewed, mode)
        kept_ensps.extend(kept)
        # count the rows removed for this gene
        if kept:
            # gene kept: its other rows (other ENSPs) are dropped
            n_extra_rows += len(ensps) - len(kept)
            if len(ensps) > len(kept):
                n_extra_genes += 1
        else:
            # gene removed: all its rows are dropped
            removed_rows_by_case[case] = removed_rows_by_case.get(case, 0) + len(ensps)
        
        for row in kept:
            if ensps.loc[row, "_status"] in ("unassigned","alternative"):
                kept_non_canonical.append(gene_name)
        
        # if the gene is removed, it is stored and grouped by the reason it was originally removed (for the summary)
        if not kept:
            if case not in removed_genes_by_case:
                removed_genes_by_case[case] = []
            removed_genes_by_case[case].append(gene_name)
    
    # The filtered table contains the kept rows, in their original order, without the helper columns and with a new column for the UniProt ID
    out = (input_df.loc[sorted(kept_ensps)]
           .drop(columns=["_status", "_n_values"])
           .rename(columns={"_accession": "UniProt_ID"}))
    
    # summary of the filtering, in the order of the steps
    print(f"\nUniProt isoform mapping summary ({mode} mode)")
    print(f"  Input ENSP IDs (rows):  {n_input_rows} rows")
    print(f"  Input genes:  {n_input_genes} genes")
    print(f"  - {n_no_name} ENSP IDs removed: no gene name was found in the dataset")
    print(f"  - {len(gene_name_mismatch)} ENSP IDs removed ({n_mismatch_genes} genes affected): "
          f"UniProt assigns the ENSP to a different gene than the dataset (see warnings above)")
    for case, gene_names in removed_genes_by_case.items():
        print(f"  - {removed_rows_by_case[case]} ENSP IDs removed ({len(gene_names)} genes removed): "
              f"{CASE_DESCRIPTIONS[case]}")
    
    print(f"  - {n_extra_rows} ENSP IDs removed ({n_extra_genes} genes affected): "
          f"number of ENSPs removed for the kept genes (since only one ensp is kept per protein)")
    print(f"  Output: {len(out)} rows; {out['Name'].nunique()} genes; {out['UniProt_ID'].nunique()} proteins")
    # check that every removed row is accounted for
    n_removed = n_no_name + len(gene_name_mismatch) + sum(removed_rows_by_case.values()) + n_extra_rows
    if n_input_rows - n_removed != len(out):
        print("  Warning: the removed rows do not add up to the difference between input and output")

    if kept_non_canonical:
        print(f"  Note: {len(kept_non_canonical)} of the kept genes have no confirmed canonical ENSP "
              f"(ENSP linked to unassigned or alternative isoform only). "
              f"Examples: {', '.join(kept_non_canonical[:10])}\n")
    
    # Put the original column names back and restore the original layout
    out = out.reset_index(drop=True)
    out = out.rename(columns={"Name": gene_name_col, "Database_ID": database_id_col})
    if index_names:
        out = out.set_index(index_names + ["UniProt_ID"])

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
     # transpose the abundance table to get proteins as rows and samples as columns
    try:
        proteomic_dataset_transpose = proteomic_cancer_type_cptac.T
        print(f"Total number of entries in the dataset: {len(proteomic_dataset_transpose)}")
        # canonical isoform mapping 
        matrix_canonical_isoform = map_uniprot_isoforms(proteomic_dataset_transpose, mode=args.isoform_mode, use_version=args.version_id) 
        # drop rows (proteins) with at least one NaN value across samples
        proteomic_dataset_transpose_clean = matrix_canonical_isoform.dropna(axis =0, how = 'any')


        # report how many rows and genes were removed because of missing values
        n_rows_before = len(matrix_canonical_isoform)
        n_rows_after = len(proteomic_dataset_transpose_clean)
        n_genes_before = matrix_canonical_isoform.index.get_level_values(0).nunique()
        n_genes_after = proteomic_dataset_transpose_clean.index.get_level_values(0).nunique()
        print(f"Missing values filter: {n_rows_before - n_rows_after} rows "
              f"({n_genes_before - n_genes_after} genes) removed because they have at least one missing value")
        print(f"Final matrix: {n_rows_after} rows ({n_genes_after} genes)")
        
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