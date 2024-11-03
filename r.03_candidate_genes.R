library(ggvenn)

setwd('/data/home/jiangminghui/project-0079/03_candidate_genes/')

res <-readRDS('../01_DEGs/01_DEGs_res.rds')

res$significant <- ifelse(res$adj.P.Val < 0.05 & res$logFC > 0.5, "Upregulated",
                          ifelse(res$adj.P.Val < 0.05 & res$logFC < -0.5, "Downregulated", "Not Significant"))

a <-readRDS('../00_bulk_rawdata/Macro_genes.rds')

b <- subset(res,significant !="Not Significant") %>%rownames()


intersect(a,b)
Overlapgenes <-intersect(a,b)
saveRDS(Overlapgenes,'./02_Overlapgenes.rds')
saveRDS(Overlapgenes,'./Overlapgenes.rds')


pdf('./venn_diagram.pdf',width = 5,height = 3)
ggvenn(list(Macro_genes = a, Significant_genes = b),
       fill_color = c("#F2A7AF", "#C7621A"))
dev.off()

png('./venn_diagram.png',width = 6,height = 4,res = 300,units = 'in')
ggvenn(list(Macro_genes = a, Significant_genes = b),
       fill_color = c("#F2A7AF", "#C7621A"))
dev.off()



length(Overlapgenes)
DEGs_entrez <- mapIds(org.Hs.eg.db, keys = Overlapgenes, column = "ENTREZID", keytype = "SYMBOL", multiVals = "first")
head(DEGs_entrez)

# ego <- enrichGO(gene          = DEGs_entrez,
#                 #universe      = names(geneList),
#                 OrgDb         = org.Hs.eg.db,
#                 ont           = "ALL",
#                 pAdjustMethod = "BH",
#                 pvalueCutoff  = 0.01,
#                 qvalueCutoff  = 0.01,
#                 readable      = TRUE)
# 
# barplot(ego, showCategory=10, title="GO Enrichment Bar Plot")
# dotplot(ego, showCategory=10, title="GO Enrichment Dot Plot")
#         
# pdf(file = "01_GO_Enrichment_Barplot.pdf")
# barplot(ego, showCategory=10, title="GO Enrichment Bar Plot")
# dev.off()
# 
# png(file = "01_GO_Enrichment_Barplot.png", width = 8, height = 6, units = 'in', res = 300)
# barplot(ego, showCategory=10, title="GO Enrichment Bar Plot")
# dev.off()
# 
# # 保存点图
# pdf(file = "01_GO_Enrichment_Dotplot.pdf")
# dotplot(ego, showCategory=10, title="GO Enrichment Dot Plot")
# dev.off()
# 
# png(file = "01_GO_Enrichment_Dotplot.png", width = 8, height = 6, units = 'in', res = 300)
# dotplot(ego, showCategory=10, title="GO Enrichment Dot Plot")
# dev.off()


ego_BP <- enrichGO(gene = DEGs_entrez, OrgDb = org.Hs.eg.db, ont = "BP", pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.05, readable = TRUE)
ego_CC <- enrichGO(gene = DEGs_entrez, OrgDb = org.Hs.eg.db, ont = "CC", pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.05, readable = TRUE)
ego_MF <- enrichGO(gene = DEGs_entrez, OrgDb = org.Hs.eg.db, ont = "MF", pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.05, readable = TRUE)

# 绘制并保存BP的柱状图
pdf(file = "01_BP_Top5_Barplot.pdf", width = 8, height = 4)
barplot(ego_BP, showCategory = 5, title = "Top 5 BP Enrichment Bar Plot")
dev.off()

png(file = "01_BP_Top5_Barplot.png", width = 8, height = 4, units = 'in', res = 300)
barplot(ego_BP, showCategory = 5, title = "Top 5 BP Enrichment Bar Plot")
dev.off()

# 绘制点图
pdf(file = "01_BP_Top5_Dotplot.pdf", width = 8, height = 4)
dotplot(ego_BP, showCategory = 5, title = "Top 5 BP Enrichment Dot Plot")
dev.off()

png(file = "01_BP_Top5_Dotplot.png", width = 8, height = 4, units = 'in', res = 300)
dotplot(ego_BP, showCategory = 5, title = "Top 5 BP Enrichment Dot Plot")
dev.off()


pdf(file = "01_CC_Top5_Barplot.pdf", width = 8, height = 4)
barplot(ego_CC, showCategory = 5, title = "Top 5 CC Enrichment Bar Plot")
dev.off()

png(file = "01_CC_Top5_Barplot.png", width = 8, height = 4, units = 'in', res = 300)
barplot(ego_CC, showCategory = 5, title = "Top 5 CC Enrichment Bar Plot")
dev.off()

# 绘制点图
pdf(file = "01_CC_Top5_Dotplot.pdf", width = 8, height = 4)
dotplot(ego_CC, showCategory = 5, title = "Top 5 CC Enrichment Dot Plot")
dev.off()

png(file = "01_CC_Top5_Dotplot.png", width = 8, height = 4, units = 'in', res = 300)
dotplot(ego_CC, showCategory = 5, title = "Top 5 CC Enrichment Dot Plot")
dev.off()


# 绘制柱状图
pdf(file = "01_MF_Top5_Barplot.pdf", width = 8, height = 4)
barplot(ego_MF, showCategory = 5, title = "Top 5 MF Enrichment Bar Plot")
dev.off()

png(file = "01_MF_Top5_Barplot.png", width = 8, height = 4, units = 'in', res = 300)
barplot(ego_MF, showCategory = 5, title = "Top 5 MF Enrichment Bar Plot")
dev.off()

# 绘制点图
pdf(file = "01_MF_Top5_Dotplot.pdf", width = 8, height = 4)
dotplot(ego_MF, showCategory = 5, title = "Top 5 MF Enrichment Dot Plot")
dev.off()

png(file = "01_MF_Top5_Dotplot.png", width = 8, height = 4, units = 'in', res = 300)
dotplot(ego_MF, showCategory = 5, title = "Top 5 MF Enrichment Dot Plot")
dev.off()




kegg_ego <- enrichKEGG(gene         = DEGs_entrez,
                       organism     = 'hsa',
                       pvalueCutoff = 0.05)

# 绘制柱状图
barplot(kegg_ego, showCategory = 10, title = "KEGG Enrichment Bar Plot")

# 绘制点图
dotplot(kegg_ego, showCategory = 10, title = "KEGG Enrichment Dot Plot")

# 保存柱状图和点图为PDF和PNG
# 保存柱状图
pdf(file = "02_KEGG_Enrichment_Barplot.pdf", width = 8, height = 7)
barplot(kegg_ego, showCategory = 10, title = "KEGG Enrichment Bar Plot")
dev.off()

png(file = "02_KEGG_Enrichment_Barplot.png", width = 8, height = 7, units = 'in', res = 300)
barplot(kegg_ego, showCategory = 10, title = "KEGG Enrichment Bar Plot")
dev.off()

# 保存点图
pdf(file = "02_KEGG_Enrichment_Dotplot.pdf")
dotplot(kegg_ego, showCategory = 10, title = "KEGG Enrichment Dot Plot")
dev.off()

png(file = "02_KEGG_Enrichment_Dotplot.png", width = 8, height = 6, units = 'in', res = 300)
dotplot(kegg_ego, showCategory = 10, title = "KEGG Enrichment Dot Plot")
dev.off()


top5_BP <- head(ego_BP@result$Description, 5)
top5_CC <- head(ego_CC@result$Description, 5)
top5_MF <- head(ego_MF@result$Description, 5)
top5_KEGG <- head(kegg_ego@result$Description, 10)

# 打印结果
print("前5个BP显著通路:")
print(top5_BP)

print("前5个CC显著通路:")
print(top5_CC)

print("前5个MF显著通路:")
print(top5_MF)

print("前5个KEGG显著通路:")
print(top5_KEGG)



bp_count <- nrow(ego_BP@result)
cc_count <- nrow(ego_CC@result)
mf_count <- nrow(ego_MF@result)
