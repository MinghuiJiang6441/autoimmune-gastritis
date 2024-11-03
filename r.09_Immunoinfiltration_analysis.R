rm(list = ls())
gc()
setwd('/data/home/jiangminghui/project-0079/09_Immunoinfiltration_analysis/')
expr <-readRDS('../00_bulk_rawdata/GSE233973_exp.rds') %>%log2()
group <-readRDS('../00_bulk_rawdata/GSE233973_group.rds')
library(xCell)
library(ggplot2)
library(reshape2)
library(plyr)
library(dplyr)
library(ggpubr)
library(corrplot)
xcell_results <- xCellAnalysis(expr)

# Add group information to xCell results
xcell_df <- as.data.frame(xcell_results)
xcell_df$CellType <- rownames(xcell_df)
# Melt the data for plotting
xcell_melted <- melt(xcell_df, id.vars = "CellType", variable.name = "Sample", value.name = "Score")

# Add group information to the melted data
sample_group <- data.frame(Sample = colnames(xcell_results), Group = group)
xcell_melted <- merge(xcell_melted, sample_group, by = "Sample")

# Normalize the scores to sum to 1 for each sample
xcell_melted <- ddply(xcell_melted, .(Sample), transform, Score = Score / sum(Score))
#xcell_df$Group <- group

# Melt the data for statistical testing and plotting
#xcell_melted <- melt(xcell_df, id.vars = "Group", variable.name = "CellType", value.name = "Score")
p1 <-ggplot(xcell_melted, aes(x = Sample, y = Score, fill = CellType)) +
  geom_bar(stat = "identity") +
  facet_wrap(~ Group, scales = "free_x") +
  theme(axis.text.x = element_blank(), 
        axis.ticks.x = element_blank(), 
        axis.line.x = element_blank()) +
  labs(title = "Immune Cell Infiltration in Disease vs. Normal", x = "Sample", y = "Proportion", fill = "Cell Type") +
  theme_minimal() +
  theme(legend.position = "bottom",
        legend.text = element_text(size = 6),
        legend.title = element_text(size = 8),
        legend.key.size = unit(0.5, "cm"))
print(p1)        
ggsave(filename = './01_immunoinfiltration.pdf',plot = p1 ,width = 17,height = 10)

ggsave(filename = './01_immunoinfiltration.png',plot = p1 ,width = 17,height = 10)

xcell_df <- as.data.frame(t(xcell_results))
xcell_df$Group <- group
xcell_df$Sample <- rownames(xcell_df)
# xcell_df$Group <- rep(group, each = 1)

# Melt the data for statistical testing and plotting
xcell_melted <- melt(xcell_df, id.vars = c("Sample", "Group"), variable.name = "CellType", value.name = "Score")
dim(xcell_df)
#
p_values <- xcell_melted %>%  
  group_by(CellType) %>%  
  dplyr::summarise(p_value = wilcox.test(Score ~ Group)$p.value)

# Filter for significantly different cell types (p < 0.05)
significant_cells <- p_values %>%
  filter(p_value < 0.05) %>%
  pull(CellType)
saveRDS(significant_cells,'./02_significant_cells.rds')

p_values_df <- data.frame(CellType = p_values$CellType, p_value = p_values$p_value)

write.csv(p_values_df, "02_cell_types_p_values.csv", row.names = FALSE)

p2 <-ggplot(xcell_melted %>% filter(CellType %in% significant_cells), aes(x = CellType, y = Score, fill = Group)) +
  geom_boxplot(width = 0.5, position = position_dodge(0.9), size = 0.2, outlier.size = 0.5, outlier.shape = NA) +  # Exclude outliers
  labs(title = "Significant Immune Cell Infiltration Differences (p < 0.05)", x = "Cell Type", y = "Score") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  stat_compare_means(aes(group = Group), label = "p.signif", method = "wilcox.test") +
  scale_fill_manual(values = c("#F20F22", "#0455BF")) +
  ylim(0, 0.24)  

ggsave('./02_Significant_Immune_Cell.pdf',plot = p2,width = 14,height = 8)

ggsave('./02_Significant_Immune_Cell.png',plot = p2,width = 14,height = 8)


xcell_df <- as.data.frame(t(xcell_results))
xcell_significant <- xcell_df[, significant_cells]

cor_matrix <- cor(xcell_significant, method = "spearman")

# Perform correlation test for significant immune cells
cor_test <- cor.mtest(xcell_significant, conf.level = 0.95)

# Filter for correlations with r > 0.3 and p < 0.05 for significant immune cells
cor_filtered <- cor_matrix
cor_filtered[abs(cor_matrix) <= 0.3 | cor_test$p > 0.05] <- NA

# Define custom colors for the correlation heatmap
col <- colorRampPalette(c("#0455BF", "white", "#F20F22"))(200)

pdf('./03_Spearman_Correlation_Significant_Immune_Cells.pdf',width =8,height =8)
corrplot(cor_filtered, method = "color", type = "upper", tl.col = "black", na.label = " ",
         col = col, tl.cex = 0.7, tl.srt = 45)
dev.off()

png('./03_Spearman_Correlation_Significant_Immune_Cells.png',width =8,height =8,units = 'in',res = 300)
corrplot(cor_filtered, method = "color", type = "upper", tl.col = "black", na.label = " ",
         col = col, tl.cex = 0.7, tl.srt = 45)
dev.off()

max_cor <- which(cor_matrix == max(cor_matrix[upper.tri(cor_matrix, diag = FALSE)]), arr.ind = TRUE)
cell_type_1 <- rownames(cor_matrix)[max_cor[1]]
cell_type_2 <- colnames(cor_matrix)[max_cor[2]]

cat(cell_type_1,cell_type_2)

highest_correlation <- cor_matrix[max_cor]
highest_p_value <- cor_test$p[max_cor]
# 获取当前目录的路径
current_dir <- getwd()

# 获取当前目录的文件夹名字
dir_name <- basename(current_dir)

# 保存当前环境到文件
save.image(file = paste0(current_dir, "/", 'Rdata_',dir_name, ".RData"))
