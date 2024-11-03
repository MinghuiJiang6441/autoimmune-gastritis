rm(list = ls())
gc()
library(GSVA)
library(GSEABase)
library(ggcor)
library(ggplot2)
library(ggtext)
library(xCell)

setwd('/data/home/jiangminghui/project-0079/10_Biologicalprocess_regulation/')
expr <-readRDS('../00_bulk_rawdata/GSE233973_exp.rds') %>%log2()
group <-readRDS('../00_bulk_rawdata/GSE233973_group.rds')
feature_gene <-readRDS('../05_machine_learning/03_feature_gene.rds')

gmt_file <- './c5.all.v7.0.symbols.gmt'
gmt_data <- getGmt(gmt_file)
str(gmt_data)
pathway_names <- names(gmt_data)
macrophage_related <- pathway_names[grepl("macrophage", pathway_names, ignore.case = TRUE)]
macrophage_pathways <- gmt_data[macrophage_related]
macrophage_pathways_list <- lapply(macrophage_pathways, geneIds)
names(macrophage_pathways_list) <-macrophage_related
xcell_results <- xCellAnalysis(expr)
significant_cells <-readRDS('../09_Immunoinfiltration_analysis/02_significant_cells.rds')
xcell_df <- as.data.frame(t(xcell_results))
xcell_significant <- xcell_df[, significant_cells]

expr <-as.matrix(expr)
gsva_scores <- gsva(expr, macrophage_pathways_list, method = "gsva", kcdf = "Gaussian")
# dim(expr)
# dim(gsva_scores)
# gsva_scores <-rbind(expr[feature_gene[1], ],gsva_scores)
# rownames(gsva_scores)[1] <-feature_gene[1]

# cor_results <- cor(t(gsva_scores), method = "spearman")

expr1 <-expr[feature_gene[1],] %>% as.data.frame()

dim(expr1)
dim(gsva_scores)
mantel <- mantel_test( expr1,t(gsva_scores), mantel.fun = 'mantel.randtest',spec.dist.method = 'bray', env.dist.method = 'euclidean'
                      ) %>% 
  mutate(r_value = cut(r, breaks = c(-Inf, 0.25, 0.5, Inf), 
                       labels = c('<0.25', '0.25-0.5', '>=0.5'), right = FALSE),
         p_value = cut(p.value, breaks = c(-Inf, 0.001, 0.01, 0.05, Inf), 
                       labels = c('<0.001', '0.001-0.01', '0.01-0.05', '>=0.05'), right = FALSE))

gsva_scores1 <-rbind(expr[feature_gene[1], ],gsva_scores)
rownames(gsva_scores)[1] <-feature_gene[1]
cor_results <- cor(t(gsva_scores), method = "spearman")
mantel$spec <-rep(feature_gene[1])
# 
# p <- quickcor(cor_results, type = "lower") +
#   geom_square() +
#   anno_link(aes(colour = p_value, size = r_value), data = mantel) +
#   scale_size_manual(values = c(0.5, 1, 2)) +
#   scale_colour_manual(values = c("#D95F02", "#1B9E77", "#A2A2A288")) +
#   scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0, space = "Lab", na.value = "grey50") +  # 重新定义颜色
#   guides(size = guide_legend(title = "Mantel's r",
#                              override.aes = list(colour = "grey35"),
#                              order = 2),
#          colour = guide_legend(title = "Mantel's p",
#                                override.aes = list(size = 3),
#                                order = 1),
#          fill = guide_colorbar(title = "Pearson's r", order = 3)) +
#   theme(axis.text.x = element_blank() )+#添加图例
#   theme(legend.position = "right",
#         axis.text.y = element_text(size = 4, color = "black"),
#      panel.background= element_blank(),#无背景
#      axis.text= element_markdown(color="black", size= 3)
#      )
# 
# 
# 
# mantel <- mantel_test( expr1,xcell_significant, mantel.fun = 'mantel.randtest',spec.dist.method = 'bray', env.dist.method = 'euclidean'
# ) %>%
#   mutate(r_value = cut(r, breaks = c(-Inf, 0.25, 0.5, Inf),
#                        labels = c('<0.25', '0.25-0.5', '>=0.5'), right = FALSE),
#          p_value = cut(p.value, breaks = c(-Inf, 0.001, 0.01, 0.05, Inf),
#                        labels = c('<0.001', '0.001-0.01', '0.01-0.05', '>=0.05'), right = FALSE))
# xcell_significant1 <-rbind(expr[feature_gene[1], ],xcell_significant)
# rownames(xcell_significant1)[1] <-feature_gene[1]
# cor_results <- cor(xcell_significant1, method = "spearman")
# mantel$spec <-rep(feature_gene[1])
# quickcor(cor_results, type = "upper") +
#   geom_square() +
#   anno_link(aes(colour = p_value, size = r_value), data = mantel) +
#   scale_size_manual(values = c(0.5, 1, 2)) +
#   scale_colour_manual(values = c("#D95F02", "#1B9E77", "#A2A2A288")) +
#   scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0, space = "Lab", na.value = "grey50") +  # 重新定义颜色
#   guides(size = guide_legend(title = "Mantel's r",
#                              override.aes = list(colour = "grey35"),
#                              order = 2),
#          colour = guide_legend(title = "Mantel's p",
#                                override.aes = list(size = 3),
#                                order = 1),
#          fill = guide_colorbar(title = "Pearson's r", order = 3)) +
#   theme(axis.text.x = element_blank(),)+#添加图例
#   theme(
#     legend.position = "right",
#     axis.text.y = element_text(size = 6, color = "black"),
#     panel.background= element_blank(),#无背景
#     axis.text= element_markdown(color="black", size= 5)
#   )


gsva_scores <- gsva(expr, macrophage_pathways_list, method = "gsva", kcdf = "Gaussian")

generate_plots <- function(expr, feature_gene, macrophage_pathways_list, xcell_significant, output_prefix) {
  # 计算GSVA分数
 # gsva_scores <- gsva(expr, macrophage_pathways_list, method = "gsva", kcdf = "Gaussian")
  
  # 循环处理每个基因
  for (gene in feature_gene) {
    expr1 <- expr[gene, ] %>% as.data.frame()
    mantel <- mantel_test(expr1, t(gsva_scores), mantel.fun = 'mantel.randtest', 
                          spec.dist.method = 'bray', env.dist.method = 'euclidean') %>%
      mutate(r_value = cut(r, breaks = c(-Inf, 0.25, 0.5, Inf),
                           labels = c('<0.25', '0.25-0.5', '>=0.5'), right = FALSE),
             p_value = cut(p.value, breaks = c(-Inf, 0.001, 0.01, 0.05, Inf),
                           labels = c('<0.001', '0.001-0.01', '0.01-0.05', '>=0.05'), right = FALSE))
    
    gsva_scores1 <- rbind(expr[gene, ], gsva_scores)
    rownames(gsva_scores1)[1] <- gene
    cor_results <- cor(t(gsva_scores1), method = "spearman")
    mantel$spec <- rep(gene, length(mantel$p_value))
    
    # 绘制第一个图
    p1 <- quickcor(cor_results, type = "lower") +
      geom_square() +
      anno_link(aes(colour = p_value, size = r_value), data = mantel) +
      scale_size_manual(values = c(0.5, 1, 2)) +
      scale_colour_manual(values = c("#D95F02", "#1B9E77", "#A2A2A288")) +
      scale_fill_gradient2(low = "#6E7348", mid = "white", high = "#8C3D20", midpoint = 0, space = "Lab", na.value = "grey50") +
      guides(size = guide_legend(title = "Mantel's r",
                                 override.aes = list(colour = "grey35"),
                                 order = 2),
             colour = guide_legend(title = "Mantel's p",
                                   override.aes = list(size = 3),
                                   order = 1),
             fill = guide_colorbar(title = "Pearson's r", order = 3)) +
      theme(axis.text.x = element_blank() )+#添加图例
      theme(legend.position = "right",
            axis.text.y = element_text(size = 4, color = "black"),
            panel.background= element_blank(),#无背景
            axis.text= element_markdown(color="black", size= 3)
      )
    
    
    # 保存第一个图
    ggsave(filename = paste0(output_prefix, "_", gene, "_lower.pdf"), plot = p1,width = 7,height = 9)
    ggsave(filename = paste0(output_prefix, "_", gene, "_lower.png"), plot = p1,width = 7,height = 9)
    
    # 处理xcell_significant数据
    mantel <- mantel_test(expr1, xcell_significant, mantel.fun = 'mantel.randtest',
                          spec.dist.method = 'bray', env.dist.method = 'euclidean') %>%
      mutate(r_value = cut(r, breaks = c(-Inf, 0.25, 0.5, Inf),
                           labels = c('<0.25', '0.25-0.5', '>=0.5'), right = FALSE),
             p_value = cut(p.value, breaks = c(-Inf, 0.001, 0.01, 0.05, Inf),
                           labels = c('<0.001', '0.001-0.01', '0.01-0.05', '>=0.05'), right = FALSE))
    
    xcell_significant1 <- rbind(expr[gene, ], xcell_significant)
    rownames(xcell_significant1)[1] <- gene
    cor_results <- cor(xcell_significant1, method = "spearman")
    mantel$spec <- rep(gene, length(mantel$p_value))
    
    # 绘制第二个图
    p2 <- quickcor(cor_results, type = "upper") +
      geom_square() +
      anno_link(aes(colour = p_value, size = r_value), data = mantel) +
      scale_size_manual(values = c(0.5, 1, 2)) +
      scale_colour_manual(values = c("#D95F02", "#1B9E77", "#A2A2A288")) +
      scale_fill_gradient2(low = "#6E7348", mid = "white", high = "#8C3D20", midpoint = 0, space = "Lab", na.value = "grey50") +
      guides(size = guide_legend(title = "Mantel's r",
                                 override.aes = list(colour = "grey35"),
                                 order = 2),
             colour = guide_legend(title = "Mantel's p",
                                   override.aes = list(size = 3),
                                   order = 1),
             fill = guide_colorbar(title = "Pearson's r", order = 3)) +
      theme(axis.text.x = element_blank(),)+#添加图例
      theme(
        legend.position = "right",
        axis.text.y = element_text(size = 6, color = "black"),
        panel.background= element_blank(),#无背景
        axis.text= element_markdown(color="black", size= 5)
      )
    
    # 保存第二个图
    ggsave(filename = paste0(output_prefix, gene, "_upper.pdf"), plot = p2,width = 7,height = 7)
    ggsave(filename = paste0(output_prefix, gene,"_upper.png"), plot = p2,width = 7,height = 7)
  }
}
  

   generate_plots(expr, feature_gene, macrophage_pathways_list, xcell_significant, "output_plot")

   
   cat(feature_gene,sep = '\n')
   
   
   
   
   
   