rm(list = ls())
gc()
setwd('/data/home/jiangminghui/project-0079/08_GSEAanalysis/')
expr <-readRDS('../00_bulk_rawdata/GSE233973_exp.rds') %>% log2()
feature_gene <-readRDS('../05_machine_learning/03_feature_gene.rds')
.libPaths("/data/nas2/Software/miniconda3/envs/public_R/lib/R/library"  )

library(limma)
library(clusterProfiler)
library(GSEABase)
library(GSVA)
library(msigdbr)
library(enrichplot)

msigdb <- msigdbr(species = "Homo sapiens", category = "C2")
kegg_gene_sets <- msigdb %>% dplyr::select(gs_name, gene_symbol)
expr <- as.matrix(expr)


perform_gsea <- function(gene) {
  # Calculate median value
  median_value <- expr[gene, ] %>% as.numeric() %>% median()
  
  # Divide samples into high and low expression groups
  high_expression_group <- expr[,which(expr[gene, ] > median_value)]
  low_expression_group <- expr[,which(expr[gene, ] <= median_value)]
  
  # Perform differential expression analysis using limma
  design <- model.matrix(~ 0 + factor(c(rep("High", ncol(high_expression_group)), rep("Low", ncol(low_expression_group)))))
  colnames(design) <- c("High", "Low")
  fit <- lmFit(cbind(high_expression_group, low_expression_group), design)
  contrast_matrix <- makeContrasts(High - Low, levels=design)
  fit <- contrasts.fit(fit, contrast_matrix)
  fit <- eBayes(fit)
  top_table <- topTable(fit, adjust="fdr", number=nrow(fit))
  
  # Filter and sort genes
  filtered_genes <- top_table[top_table$adj.P.Val < 0.05, ]
  sorted_genes <- filtered_genes[order(-filtered_genes$logFC), ]$logFC
  names(sorted_genes) <- rownames(filtered_genes[order(-filtered_genes$logFC), ])
  
  # Perform GSEA using clusterProfiler
  gsea_results <- GSEA(sorted_genes, TERM2GENE=kegg_gene_sets, pvalueCutoff=0.05, eps = 0)
  if (is.null(gsea_results) || nrow(gsea_results) == 0) {
    return(NULL)
  }
  significant_pathways <- gsea_results[gsea_results$NES > 1, ]
  top5_pathways <- significant_pathways[order(significant_pathways$p.adjust), ][1:5, ]
  # Plot the top 5 pathways
  pathway <- top5_pathways[1:5, "ID"]
  plot <- gseaplot2(gsea_results, geneSetID=pathway, title=paste(gene, "Top 5 Pathways"))
  
  # Save plots as PDF and PNG
  pdf(file = paste0("GSEA_", gene, ".pdf"))
  print(plot)
  dev.off()
  
  png(file = paste0("GSEA_", gene, ".png"), width=8, height=6,units = 'in',res = 300)
  print(plot)
  dev.off()
  
  return(plot)
  return(All_Top5_Pathways)
}

plot_list <- list()

for (gene in feature_gene) {
  plot_list[[gene]] <- perform_gsea(gene)
  
}
