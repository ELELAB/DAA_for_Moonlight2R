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
  
  limma_results_joined <- limma_results

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
  make_option(c("-p", "--plot"), action = "store_true", default = FALSE, help="If the flag is specified, plots will be generated"),
  make_option(c("-t", "--threshold_pval"), type= "numeric", default = 0.05, help="If the flag is specified, a custom p value threshold will be applied to the volcano plot")
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

if(any(duplicated(data_matrix$Name)))
{
  stop("The input matrix contains duplicate gene names.")
}

# run limma
limma_input_matrix <- prepareLimmaInput(data_matrix)
limma_results <- runLimma(limma_input_matrix, sample_labels_df, contrast_labels)

# remove Avg expression, t-statistics and Bonferroni correction columns from the output df 
limma_results <- limma_results[, !names(limma_results) %in% c("AveExpr","t","B")]
# assign gene names as rownames
limma_results <- column_to_rownames(limma_results, var = "Name")
limma_results <- limma_results %>% rename(PValue = P.Value)

cancer_type <- args$cancer_type

results_dir_path <- file.path("limma_results")
if(!dir.exists("limma_results"))
{
  dir.create("limma_results", recursive = TRUE)
}

if(cancer_type != "None")
{
  write.csv(limma_results, paste0(results_dir_path,"/limma_table_", cancer_type,".csv"), row.names = TRUE)
} else
{
  write.csv(limma_results, paste0(results_dir_path,"/limma_table.csv"), row.names = TRUE)
}
  
if(args$plot)
  plotLimmaResults(limma_results, args$threshold_pval, cancer_type)

