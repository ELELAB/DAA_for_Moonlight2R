library(limma)
library(optparse)
library(dplyr)
library(tidyverse)
#library(ggplot2)

# --- FUNCTIONS ---
# load abundance matrix file
load_data <- function(data_path) {
  original_data <- read.csv(data_path, check.names = FALSE)
  return (original_data)
}

mapProteinIsoforms <- function(abundance_matrix) {
  map_df <- abundance_matrix[, c("Name","Database_ID")]
  return(map_df)
}

prepareLimmaInput <- function(abundance_matrix) {
  limma_input_matrix <- abundance_matrix[, !names(abundance_matrix) %in% c("Name")] # keep only database IDs (keep isoform level)
  rownames(limma_input_matrix) <- limma_input_matrix[[1]]
  limma_input_matrix[[1]] <- NULL
  return(limma_input_matrix)
}

groupSamples <- function(limma_input_matrix) {
  
  sample_IDs <- colnames(limma_input_matrix)
  df_out <- data.frame(Sample_IDs = sample_IDs, Group_Labels = NA)
  for (s in sample_IDs)
  {
    is_normal_sample <- grepl("\\.N$", as.character(s))
    #print(s)
    if (is_normal_sample)
    {
      df_out$Group_Labels[df_out$Sample_IDs == s]<- "N"
    }
    else 
    {
      df_out$Group_Labels[df_out$Sample_IDs == s]<- "T"
    }
    all_same <- n_distinct(df_out$Group_Labels) == 1
    if (all_same)
    {
      stop("At least two samples group should be defined. Only one found. Exiting")
    }
  }
  return (df_out)
}

runLimma <- function(matrix, groups) {

  groups <- factor(groups$Group_Labels)
  design <- model.matrix(~0 + groups)

  fit <- lmFit(matrix, design)
  
  design_cols <- colnames(design)
  contrast_exp <- paste0(design_cols[2], " - ", design_cols[1])
  contrasts_matrix <- makeContrasts(contrasts=contrast_exp, levels = design)

  fit2 <- contrasts.fit(fit, contrasts_matrix)

  fit3 <- eBayes(fit2, robust = TRUE)

  limma_results <- topTable(fit3, coef = 1, adjust.method="BH", number=Inf) %>% 
                  rownames_to_column("Database_ID") %>% 
                  as_tibble()
  limma_results_joined <- protein_map_df %>% left_join(limma_results, by = "Database_ID")

  limma_results_joined$Direction <- ifelse(limma_results_joined$logFC > 0, "up", "down") %>% as.factor()

  return(limma_results_joined)

}

plotLimmaResults <- function(limma_table, threshold_pval, cancer_type) {

  if(cancer_type != "None"){
    filename <- paste0("limma_volcano_plot_", cancer_type, ".png")
    png(filename)
  } else {
    png("limma_volcano_plot.png")
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

# --- RUN LIMMA --- 
# parse argmuments
option_list <- list(
  make_option(c("-m", "--matrix")),
  make_option(c("-c", "--cancer_type"), type = "character", default="None"),
  make_option(c("-p", "--plot"), action = "store_true", default = FALSE),
  make_option(c("-t", "--threshold_pval"), type="double", default=0.05)
)

parser <- OptionParser(option_list = option_list)
args <- parse_args(parser)

data_path <- args$matrix
data_matrix <- load_data(data_path)

protein_map_df <- mapProteinIsoforms(data_matrix)
limma_input_matrix <- prepareLimmaInput(data_matrix)
sample_group_df <- groupSamples(limma_input_matrix) # TODO: need to decide whether to do this in limma or in the dataset generation script. in any case we miss generalisation

limma_results <- runLimma(limma_input_matrix, sample_group_df)

cancer_type <- args$cancer_type
if(cancer_type != "None")
{
  write.csv(limma_results, paste0("limma_table_", cancer_type,".csv"))
} else { 
  write.csv(limma_results, "limma_table.csv") 
}

if(args$plot)
{
  plotLimmaResults(limma_results, args$threshold_pval, cancer_type)
}


# TODO:
# - improve plotting layout
# - error handling 
# - decide how to generalize the group labels assignment and the DATABASE IDs column as well (not always present) - may want to add some flag option for this
# - the scirpt currently handles the isoform separately, create an option to keep only the main isoform (by mapping with uniprot ID...)



