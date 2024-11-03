rm(list = ls())
gc()
setwd("/data/home/jiangminghui/project-0079/01_DEGs/")

library(ggplot2)
library(cowplot)
library(harmony)
library(readr) # 用于读取txt文件
#install.packages("scCustomize")

library(tidyr)
library(gplots)
library(reshape2)
library(DESeq2)
library(edgeR)
library(igraph)
library(tinyarray)
library(RColorBrewer)
library(circlize)
library(ComplexHeatmap)
library(ggrepel)
library(dplyr)

expr<-readRDS('../00_bulk_rawdata/GSE233973_exp.rds')
group<-readRDS('../00_bulk_rawdata/GSE233973_group.rds')


ex <- expr
qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC <- (qx[5] > 100) ||
  (qx[6]-qx[1] > 50 && qx[2] > 0) ||
  (qx[2] > 0 && qx[2] < 1 && qx[4] > 1 && qx[4] < 2)

if (LogC) { ex[which(ex <= 0)] <- NaN
expr <- log2(ex)
print("log2 transform finished")}else{print("log2 transform not needed")}

boxplot(expr,main="expression value",las=2)

design <- model.matrix(~ 0 + group)
colnames(design) <- c("disease", "normal")
contrast.matrix<-makeContrasts(paste0(unique(group),collapse = "-"),levels = design)
contrast.matrix ##这个矩阵声明，我们要把progres.组跟stable进行差异分析比较
fit <- lmFit(expr,design)
##step2
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)  ## default no trend !!!
##eBayes() with trend=TRUE
##step3
tempOutput = topTable(fit2 ,adjust.method = "fdr", coef=1, n=Inf)
res = na.omit(tempOutput) 
#write.csv(nrDEG2,"limma_notrend.results.csv",quote = F)
head(res)

saveRDS(res,'./01_DEGs_res.rds')

# 添加显著性标记
res$significant <- ifelse(res$adj.P.Val < 0.05 & res$logFC > 0.5, "Upregulated",
                          ifelse(res$adj.P.Val < 0.05 & res$logFC < -0.5, "Downregulated", "Not Significant"))

table(res$significant)
# 将res对象转换为数据框
res_df <- as.data.frame(res)

# 筛选出上调基因中log2FoldChange最大的前20个基因
top_upregulated <- res_df %>%
  filter(significant == 'Upregulated') %>%
  arrange(desc(logFC)) %>%
  head(10)

# 筛选出下调基因中log2FoldChange最小的前20个基因
top_downregulated <- res_df %>%
  filter(significant == 'Downregulated') %>%
  arrange(logFC) %>%
  head(10)

# 合并上调和下调的基因
top_genes <- bind_rows(top_upregulated, top_downregulated)
top_genes$gene <- rownames(top_genes)

# 绘制火山图
volcano_plot <- ggplot(res_df, aes(x = logFC, y = -log10(adj.P.Val), color = significant)) +
  geom_point(alpha = 0.8, size = 1.75) +
  scale_color_manual(values = c("Not Significant" = "grey", "Upregulated" = "#A60A27", "Downregulated" = "#D9A60D")) +
  theme_minimal() +
  labs(title = "Volcano Plot", x = "Log2 Fold Change", y = "-Log10 Adjusted P-value") +
  geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed", color = "blue") +  # 添加log2FoldChange阈值线
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "blue") +  # 添加p值阈值线
  xlim(-21,22) +  # 设置x轴范围
  ylim(0,26) +  # 设置y轴范围
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10),
    panel.grid.major = element_blank(),  # 去除主网格线
    panel.grid.minor = element_blank(),  # 去除次网格线
    panel.background = element_blank(),  # 去除背景
    axis.line = element_line(color = "black")  # 添加轴线
  ) +
  geom_point(size = 3, shape = 1, data = top_genes) +
  ggrepel::geom_label_repel(
    aes(label = gene), 
    data = top_genes,
    color="black"
  )


# +  geom_text_repel(data = top_genes, aes(label = gene), size = 3, max.overlaps = 10, 
#                   box.padding = 0.5, point.padding = 0.5, segment.color = 'grey50', color = "black")  # 调整标签
options(ggrepel.max.overlaps = Inf)

# 显示火山图
print(volcano_plot)

ggsave('./01_DEGs_volcano_plot.png',plot =volcano_plot,width = 10,height = 7 )
ggsave('./01_DEGs_volcano_plot.pdf',plot =volcano_plot,width = 10,height = 7 )

expr_scaled <- t(scale(t(expr[top_genes$gene,])))
group_colors <- list(Group = c("normal" = "#049DBF", "disease" = "#8C030E"))

ha <- HeatmapAnnotation(df = data.frame(Group = group), col = group_colors)

# 绘制热图
heatmap <- Heatmap(expr_scaled, 
                   name = "Expression", 
                   top_annotation = ha, 
                   show_row_names = T, 
                   show_column_names = FALSE,
                   cluster_rows = F, 
                   cluster_columns = F,
                   col = colorRamp2(c(min(expr_scaled), 0, max(expr_scaled)), c("#03178C", "white", "#A52226")))

# 显示热图
draw(heatmap)

pdf('./02_DEGs_heatmap.pdf', width = 14, height = 9)
draw(heatmap)
dev.off()

png('./02_DEGs_heatmap.png', width = 14, height = 9, units = "in", res = 300)
draw(heatmap)
dev.off()






# 获取当前目录的路径
current_dir <- getwd()

# 获取当前目录的文件夹名字
dir_name <- basename(current_dir)

# 保存当前环境到文件
save.image(file = paste0(current_dir, "/", 'Rdata_',dir_name, ".RData"))





