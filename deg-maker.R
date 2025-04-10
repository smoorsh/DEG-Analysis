#!/usr/bin/env Rscript

# Load required libraries
library(DESeq2)
library(ggplot2)
library(dplyr)
library(ggrepel)

# Define helper functions FIRST
# 1. Extract expression data for samples in a comparison group
subgem <- function(gem, anot, group) {
  subanot = subset(anot, Comparison == group)
  sample_ids = subanot$Sample
  subcounts = gem[, colnames(gem) %in% sample_ids, drop=FALSE]
  return(subcounts)
}

# 2. Filter annotations for the group of interest
subanot <- function(anot, group) {
  subanot = subset(anot, Comparison == group)
  return(subanot)
}

# 3. DESeq2 analysis function
run_deseq <- function(counts, annotation) {
  # Create DESeq dataset
  dds <- DESeqDataSetFromMatrix(countData = counts,
                                colData = annotation,
                                design = ~ Group)
  
  # Filter low count genes
  dds <- dds[rowSums(counts(dds)) >= 50,]
  
  # Run DESeq
  dds <- DESeq(dds)
  
  actual_levels <- levels(dds$Group)
  reference_level <- actual_levels[1]
  treatment_level <- actual_levels[2]
  
  # Get normalized counts
  norm <- fpm(dds)
  conditionA <- colnames(norm)[dds$Group == reference_level]
  conditionB <- colnames(norm)[dds$Group == treatment_level]
  norm_subset <- norm[, c(conditionA, conditionB)]
  
  # Get results
  res <- results(dds, name=paste0("Group_", treatment_level, "_vs_", reference_level))
  
  # Add normalized expression data
  norm_means <- data.frame(
    reference_mean = rowMeans(norm_subset[, conditionA, drop=FALSE]),
    treatment_mean = rowMeans(norm_subset[, conditionB, drop=FALSE])
  )
  colnames(norm_means) <- c(paste0(reference_level, "_mean"), paste0(treatment_level, "_mean"))
  
  res_df <- as.data.frame(res)
  res_with_norm <- cbind(res_df, norm_means)
  
  return(res_with_norm)
}

# 4. Main execution function
main <- function(countfile, anotfile, outfile, add_comparison_to_filename = FALSE) {
  # Read input files
  counts = read.delim(countfile, sep=',', header=TRUE, row.names='Hugo_Symbol')
  samples = read.delim(anotfile, sep=',', row.names = NULL, check.names=FALSE)
  
  # Get unique comparison groups
  groups = unique(samples$Comparison)
  
  # Process each comparison group
  for (t in groups) {
    # Set output filename
    if (add_comparison_to_filename) {
      outname = paste0(tools::file_path_sans_ext(outfile), "_", t, ".csv")
    } else {
      outname = paste0(tools::file_path_sans_ext(outfile), ".csv")
    }
    
    # Extract data for this comparison
    subcounts = subgem(counts, samples, t)
    subannotation = subanot(samples, t)
    
    if (nrow(subannotation) >= 2) {
      results = run_deseq(subcounts, subannotation)
      f_results = subset(results, padj < 0.05)
      if (nrow(f_results) > 0) {
        o_results = f_results[order(f_results$padj),]
        write.csv(o_results, outname, row.names = TRUE)
        print(paste("Wrote", nrow(o_results), "DEGs to", outname))
      } else {
        print("No significant DEGs found for this comparison")
      }
    } else {
      print(paste("Skipping", t, "- insufficient samples"))
    }
  }
  
  # Return the output filename
  if (add_comparison_to_filename) {
    return(paste0(tools::file_path_sans_ext(outfile), "_", groups[1], ".csv"))
  } else {
    return(paste0(tools::file_path_sans_ext(outfile), ".csv"))
  }
}

# Process command line arguments
args <- commandArgs(trailingOnly = TRUE)
countfile <- args[1]
anotfile <- args[2]
outfile <- args[3]
plot_title <- args[4]
plot_output_file <- args[5]

# Run main function and get the output file path
deg_output_file <- main(countfile, anotfile, outfile, FALSE)
print(paste("Output file:", deg_output_file))

# Create and save volcano plot if the output file exists
if (file.exists(deg_output_file)) {
  deg_data <- read.csv(deg_output_file, row.names = 1)
  
  # Add significance labels
  deg_data <- deg_data %>%
    mutate(
      significance = case_when(
        padj < 0.001 & log2FoldChange > 1 ~ "Upregulated",
        padj < 0.001 & log2FoldChange < -1 ~ "Downregulated",
        TRUE ~ "Not significant"
      )
    )
  
  # Create volcano plot
  volcano_plot <- ggplot(deg_data, aes(x = log2FoldChange, y = -log10(padj))) +
    geom_point(aes(color = significance), size = 1, alpha = 0.7) +
    scale_color_manual(values = c("Upregulated" = "green", 
                                  "Downregulated" = "red", 
                                  "Not significant" = "gray")) +
    geom_hline(yintercept = -log10(0.001), linetype = "dashed", color = "gray40") +
    geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "gray40") +
    theme_minimal() +
    labs(
      title = "Volcano Plot of Differential Gene Expression",
      subtitle = plot_title,
      x = "log2 Fold Change",
      y = "-log10(adjusted p-value)",
      color = "Gene Status"
    )
  
  # Add labels for top genes
  top_genes <- deg_data %>%
    filter(padj < 0.01) %>%
    arrange(padj) %>%
    head(10)
  
  volcano_plot_with_labels <- volcano_plot +
    geom_text_repel(
      data = top_genes,
      aes(label = rownames(top_genes)),
      size = 3,
      box.padding = 0.5,
      point.padding = 0.3,
      force = 2,
      max.overlaps = 15
    )
  
  # Save plot
  ggsave(plot_output_file, volcano_plot_with_labels, width = 10, height = 8, dpi = 300)
  print(paste("Volcano plot saved as", plot_output_file))
} else {
  stop(paste("DEG output file not found:", deg_output_file))
}