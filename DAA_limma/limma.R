library(limma)
library(optparse)
library(dplyr)
library(tidyverse)
library(httr2)
library(httr)
library(jsonlite)
library(purrr)

load_data <- function(data_path) {

  # function to read a csv file- given a path

  original_data <- read.csv(data_path, check.names = FALSE)
  return (original_data)
}

getUniprotCanonicalIsoform <- function(uniprot_data){
  
  #
  # Retrieve the canonical isoform for a given Uniprot ID
  #
  
  alternative_prods <- list()
  canonical_isoform <- NA
  for (x in uniprot_data$comments)
  {
    if(x$commentType == "ALTERNATIVE PRODUCTS")
    {
      alternative_prods <- append(alternative_prods, x)
      for (prod in alternative_prods)
      {
        if(is.list(prod))
        {
          for(iso_idx in seq_along(prod))
          {
            if(is.list(prod[[iso_idx]]))
            {
              status <- prod[[iso_idx]]$isoformSequenceStatus
              if(!is.null(status))
              {
                if(prod[[iso_idx]]$isoformSequenceStatus == "Displayed")
                {
                  canonical_isoform <- prod[[iso_idx]]$isoformIds
                }
              }
            }
          }
        }
      }
    }
  }
  if (length(canonical_isoform) > 1)
    return(canonical_isoform[[1]])
  else
    return(canonical_isoform)
}

getEnsemblProtIDs <- function(uniprot_data){
  
  #
  # Retrieve the ensembl IDs for each protein isoform
  #

  ensembl_refs <- uniprot_data$uniProtKBCrossReferences
  ensembl_refs <- ensembl_refs[vapply(ensembl_refs, function(x) x$database == "Ensembl", logical(1))]
  isoform_ids <- character(length(ensembl_refs))
  ensp_ids    <- character(length(ensembl_refs))
  for(i in seq_along(ensembl_refs))
  {
    prop <- ensembl_refs[[i]]
    if(is.list(prop))
    {
      if(is.null(prop$isoformId))
      {
        next
      }
      isoform_ids[i] <- prop$isoformId
      idx <- vapply(prop$properties, function(x) x$key=="ProteinId", logical(1))
      if(any(idx))
      {
        ensp_ids[i] <- prop$properties[[which(idx)[1]]]$value
      }
      else
      {
        ensp_ids[i] <-  NA_character_
      }
    }
  }
  isoform_ensp_df <- data.frame(isoformId = isoform_ids, ensp = ensp_ids, stringsAsFactors = FALSE)
  
  return (isoform_ensp_df)
}

mapProteinIsoforms <- function(abundance_matrix) { 
  
  # Filter the input dataset to keep only Uniprot 
  # canonical isoforms It works only when an ensembl 
  # IDs column is present in the input dataset 

  char_cols <- names(abundance_matrix)[sapply(abundance_matrix, is.character)]
  if (length(char_cols) == 0)
    stop("No character columns found")
  
  # check if the abundance matrix has a Database ID column
  has_databaseID_col <- any(grepl("^database[ _]?id$", char_cols, ignore.case = TRUE))
  if(!has_databaseID_col)
  {
    print("No Database ID column found. The analysis will be run without Uniprot Isoform mapping")
    return (NULL)
  }
  map_df <- abundance_matrix[, c(char_cols[[1]], char_cols[[2]])]
  
  has_databaseID_col <- any(grepl("^database[ _]?id$", char_cols, ignore.case = TRUE))
  if(!has_databaseID_col)
  {
    print("No Database ID column found")
    return (NULL)
  }

  # get duplicate gene names- if any
  dup_names <- unique(map_df$Name[duplicated(map_df$Name)])
  # initialize a temporary df to store gene names and db ids (e.g. esembl ids)
  tmp_df <- data.frame(
  Name = character(0),
  Database_ID = character(0),
  stringsAsFactors = FALSE
  )
  # iterate over each gene entry that has duplicates
  for (gene in dup_names)
  {
    # connect to Uniprot API search to retrieve data for the current gene
    dat <- tryCatch({
      req <- request("https://rest.uniprot.org/uniprotkb/search") |>
      req_url_query(
      query = paste0("gene:", gene, " AND reviewed:true AND organism_id:9606"),
      format = "json",
      fields = "accession,gene_primary"
      )
      res <- req_perform(req)
      resp_body_json(res)

    }, 
    error = function(e) {
      message("Uniprot query failed for gene", gene, ": ", e$message)
      NULL
      }
    )
    if(is.null(dat) || length(dat$results) == 0) 
      next

    # get the uniprot ID for the current gene
    uniprot_id <- dat$results[[1]]$primaryAccession

    # retrieve uniprot data for the uniprot ID 
    dat_up <- tryCatch({

      req_up <- request(paste0("https://rest.uniprot.org/uniprotkb/", uniprot_id)) |>
      req_url_query(format = "json") |>
      req_perform()
      resp_body_json(req_up)

    },
    error = function(e){
      message("UniProt lookup failed for gene ", gene, " (UniProt ID: ", uniprot_id, "): ", e$message)
      NULL
    }
    )  
    if (is.null(dat_up)) 
      next
    # get uniprot sequence canonical isoform for the current gene
    up_canonical_isoform <- getUniprotCanonicalIsoform(dat_up)

    if (is.na(up_canonical_isoform)) 
       up_canonical_isoform <- dat_up$primaryAccession
    
    # get the ensembl protein IDs associated to each protein isoform for the current gene
    df_isoforms_ensp <- getEnsemblProtIDs(dat_up)

    if(any(df_isoforms_ensp$ensp == ""))
    {
      df_isoforms_ensp <- do.call(rbind, lapply(dat_up$uniProtKBCrossReferences, function(x) {
        if (x$database == "Ensembl" && !is.null(x$properties$proteinId)) {
          data.frame(
            isoformId = up_canonical_isoform,
            ensp = x$properties$proteinId,
            stringsAsFactors = FALSE
          )
        } else {
          NULL
        }
      }))
    }
    if(is.null(df_isoforms_ensp))
    {
      abundance_matrix <- abundance_matrix[abundance_matrix$Name != gene,]
      next
    }
    # get the ensp related to the canonical isoform only (could be more than one)
    canonical_ensp_id <- df_isoforms_ensp[df_isoforms_ensp$isoformId == up_canonical_isoform,]$ensp
    # if no ensp is associated to the uniprot canonical isoform, remove the gene from the input dataset
    if(is.null(canonical_ensp_id) || length(canonical_ensp_id) == 0)
    {
      abundance_matrix <- abundance_matrix[abundance_matrix$Name != gene,]
      next
    }
    # get the ensp ids present in the input dataframe for the current gene
    ensp_ids_map_df <- abundance_matrix$Database_ID[abundance_matrix$Name == gene]
    # filter the input dataframe to only keep the ensp ids associated to the canonical isoform
    if(length(canonical_ensp_id) > 1)
    {
      # --- if no esps associated to the canonical isoform is present in the input df- remove all the gene rows.
      if(!any(canonical_ensp_id %in% ensp_ids_map_df))
      {
        # if no ensp ids is found in the input df- this means that no canonical isoform is present for that gene and we should remove the gene entry
        abundance_matrix <- abundance_matrix[abundance_matrix$Name != gene,]
        next
      }
      else
      {
        # get only the ensp associated to the canonical isoform that are present in the input df 
        ensp_in_df <- canonical_ensp_id[canonical_ensp_id %in% map_df$Database_ID]
        # if still more than one ensp is associated to canonical isoform and present in the input df- aggregate data
        if(length(ensp_in_df) > 1)
        {
          # aggregate
          abundance_matrix <- abundance_matrix %>% filter(!(Name == gene & !Database_ID %in% ensp_in_df))
          gene_rows <- abundance_matrix %>% filter(Name == gene) %>% summarise(Name=gene,Database_ID = first(Database_ID), across(where(is.numeric), mean, na.rm = TRUE))
          abundance_matrix <- abundance_matrix %>% filter(Name != gene) %>% bind_rows(gene_rows)     
        }
        else
        {
          if(!is.na(ensp_in_df))
            abundance_matrix <- abundance_matrix[!(abundance_matrix$Name == gene & abundance_matrix$Database_ID != ensp_in_df),]
        }
      }
    }
    # if only one ensp is associated to the canonical isoform
    # keep only the row associated with it in the input df
    else if (length(canonical_ensp_id) == 1)
    {
      if(!is.na(canonical_ensp_id))
      {
        if(!canonical_ensp_id %in% abundance_matrix[abundance_matrix$Name == gene,]$Database_ID)
        {
          abundance_matrix <- abundance_matrix[abundance_matrix$Name != gene,]
          next        
        }
        else
          abundance_matrix <- abundance_matrix[!(abundance_matrix$Name == gene & abundance_matrix$Database_ID != canonical_ensp_id),]
      }
    }
  }

  return(abundance_matrix)
}

prepareLimmaInput <- function(abundance_matrix) {
  
  #
  # Manipulate the input abundance matrix to keep
  # the gene name, abundance values only, p-values and 
  #
  limma_input_matrix <- abundance_matrix[, !names(abundance_matrix) %in% c("Database_ID")]
  rownames(limma_input_matrix) <- limma_input_matrix[[1]]
  limma_input_matrix[[1]] <- NULL
  
  return(limma_input_matrix)

}


runLimma <- function(matrix, groups, contrast_labels) {
  
  #
  # Run Differential Abundance Analysis using limma-
  # given the abundance matrix and the groups/condition 
  # labels
  #

  contrast_lab_1 <- contrast_labels[1]
  contrast_lab_2 <- contrast_labels[2]
  groups <- factor(groups$Label_Group, levels = c(contrast_lab_1, contrast_lab_2))
  if(!(contrast_lab_1 %in% groups && contrast_lab_2 %in% groups ))
  {
    stop("The specified labels for the contrast matrix design do not correspond to the actual labels assigned to the samples. Please check labels definition and specify them accordingly")
  }
  design <- model.matrix(~0 + groups)
  fit <- lmFit(matrix, design)
  
  design_cols <- colnames(design)
  contrast_exp <- paste0(design_cols[1], " - ", design_cols[2])
  contrasts_matrix <- makeContrasts(contrasts=contrast_exp, levels = design)

  fit2 <- contrasts.fit(fit, contrasts_matrix)

  fit3 <- eBayes(fit2, robust = TRUE)

  limma_results <- topTable(fit3, coef = 1, adjust.method="BH", number=Inf) %>% 
                  rownames_to_column("Name") %>% 
                  as_tibble()
  
  limma_results_joined <- limma_results #protein_map_df %>% left_join(limma_results, by = "Database_ID")

  limma_results_joined$Direction <- ifelse(limma_results_joined$logFC > 0, "up", "down") %>% as.factor()

  return(limma_results_joined)

}

plotLimmaResults <- function(limma_table, threshold_pval, cancer_type) {

  #
  # Create a volcano plot from the limma output table
  #
  
  if(cancer_type != "None"){
    filename <- paste0(results_dir_path,"/limma_volcano_plot_", cancer_type, ".png")
    png(filename)
  } else {
    png(paste0(results_dir_path,"/limma_volcano_plot.png"))
  }
  limma_table$significance <- ifelse(limma_table$adj.P.Val < threshold_pval,
  "sig", "not.sig") %>%
  as.factor()
  limma_table %>% 
  ggplot(aes(x = logFC, y = -log10(adj.P.Val))) + geom_point(aes(colour = significance:Direction), size = 0.5) +
  scale_color_manual(
  values = c("black", "black", "deepskyblue", "red"), name = "",
  labels = c("Downregulated insignificant",
  "Upregulated insignificant",
  "Downregulated significant",
  "Upregulated significant")) +
  theme(axis.title.x = element_text(size = 15, vjust = -2),
  axis.title.y = element_text(size = 15, vjust = 2),
  axis.text.x = element_text(size = 12, vjust = -1),
  axis.text.y = element_text(size = 12),
  plot.background = element_rect(fill = "white"),
  panel.background = element_rect(fill = "white"),
  axis.line = element_line(linewidth = 0.5, colour = "black"),
  plot.margin = margin(10, 10, 10, 10),
  legend.position = c(0.25, 0.9)) +
  labs(x = "log2(FC)", y = "-log10(p-value)") +
  xlim(-3.1, 3.1)
}


# parse argmuments
option_list <- list(
  make_option(c("-m", "--matrix"), help="A csv file containing the protein abundance matrix"),
  make_option(c("-s", "--sample_labels", help="A csv file containing for each sample, its group/condition label")),
  make_option(c("-d", "--design_contrast"), type="character", default="T-N", help="A string indicating how to design the contrast in the differential abundance analysis 
                (e.g. if T-N is specified- the fold change for each protein will be the result of Tumor condition - Normal condition). 
                Please remember to use the same labels definition as in the input provided in --sample_labels"),
  make_option(c("-c", "--cancer_type"), type = "character", default="None", help="A string indicating the cancer type"),
  make_option(c("-p", "--plot"), action = "store_true", default = FALSE, help="If the flag is specified, plots will be generated")
)

parser <- OptionParser(option_list = option_list)
args <- parse_args(parser)

# parse data paths
data_matrix_path <- args$matrix
sample_labels_path <- args$sample_labels
# parse contrast matrix labels
contrast_labels <- strsplit(args$design_contrast, "-", fixed = TRUE)[[1]]
# load abundance matrix
data_matrix <- load_data(data_matrix_path)
# load sample labels dataframe
sample_labels_df <- load_data(sample_labels_path)

filtered_data_matrix <- mapProteinIsoforms(data_matrix)
if(is.null(filtered_data_matrix))
{
  filtered_data_matrix <- data_matrix
}
if(any(duplicated(filtered_data_matrix$Name)))
{
  stop("The input matrix contains duplicate gene names.")
}

# run limma
limma_input_matrix <- prepareLimmaInput(filtered_data_matrix)
limma_results <- runLimma(limma_input_matrix, sample_labels_df, contrast_labels)

# remove Avg expression, t-statistics and Bonferroni correction columns from the output df 
limma_results <- limma_results[, !names(limma_results) %in% c("AveExpr","t","B")]

cancer_type <- args$cancer_type

results_dir_path <- file.path("limma_results")
if(!dir.exists("limma_results"))
{
  dir.create("limma_results", recursive = TRUE)
}

if(cancer_type != "None")
{
  write.csv(limma_results, paste0(results_dir_path,"/limma_table_", cancer_type,".csv"))
} else
{
  write.csv(limma_results, paste0(results_dir_path,"/limma_table.csv"))
}
  
if(args$plot)
  plotLimmaResults(limma_results, args$threshold_pval, cancer_type)

