rm(list = ls())
gc()
setwd('/data/home/jiangminghui/project-0079/07_Correlation_functional_similarity/')
load('./Rdata_Correlation_functional_similarity.RData')
library(corrr)
library(dplyr)
library(GOSemSim)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(ggplot2)
expr <-readRDS('../00_bulk_rawdata/GSE233973_exp.rds')
feature_gene <-readRDS('../machine_learning/03_feature_gene.rds')
data <-expr[feature_gene,]
rm(expr)

x <- data %>% t()%>%
  correlate(method = 'spearman')%>% 
  rearrange() %>% 
  network_plot(min_cor = .2)

library(cor)
library(corrplot)
spearman_cor <- cor(t(data), method = "spearman", use = "complete.obs")
spearman_p <- cor.mtest(t(data), method = "spearman")

significant_cor <- spearman_cor[upper.tri(spearman_cor, diag = FALSE)]
significant_cor <- significant_cor[abs(significant_cor) > 0.3 & spearman_p$p[upper.tri(spearman_p$p, diag = FALSE)] < 0.05]

significant_cor_matrix <- matrix(NA, nrow = ncol(t(data)), ncol = ncol(t(data)))
rownames(significant_cor_matrix) <- colnames(t(data))
colnames(significant_cor_matrix) <- colnames(t(data))
significant_cor_matrix[upper.tri(significant_cor_matrix)] <- significant_cor
significant_cor_matrix[lower.tri(significant_cor_matrix)] <- t(significant_cor_matrix)[lower.tri(significant_cor_matrix)]

pdf('./01_ correlation_candidate_key_genes.pdf',width = 5,height = 5)
data %>% t()%>%
  correlate(method = 'spearman')%>% 
  rearrange() %>% 
  network_plot(min_cor = .2)
dev.off()

png('./01_ correlation_candidate_key_genes.png',width = 5,height = 5,units = 'in',res = 300)
data %>% t()%>%
  correlate(method = 'spearman')%>% 
  rearrange() %>% 
  network_plot(min_cor = .2)
dev.off()



entrez_ids <- mapIds(org.Hs.eg.db, keys = feature_gene, column = "ENTREZID", keytype = "SYMBOL") %>%as.data.frame()
entrez_ids$SYMBOL <-rownames(entrez_ids)
colnames(entrez_ids) <-c('ENTREZID','SYMBOL')
entrez_ids$ENTREZID <- as.character(entrez_ids$ENTREZID)
head(entrez_ids)

bp <- godata( 'org.Hs.eg.db', ont= "BP", computeIC = FALSE)
cc <- godata( 'org.Hs.eg.db', ont= "CC", computeIC = FALSE) 
mf <- godata( 'org.Hs.eg.db', ont= "MF", computeIC = FALSE)

simbp<- mgeneSim(entrez_ids$ENTREZID, semData= bp, measure= "Wang", drop= NULL, combine= "BMA") 
simcc<- mgeneSim(entrez_ids$ENTREZID, semData= cc, measure= "Wang", drop= NULL, combine= "BMA") 
simmf<- mgeneSim(entrez_ids$ENTREZID, semData= mf, measure= "Wang", drop= NULL, combine= "BMA")

fsim<- (simmf * simcc * simbp)^( 1/ 3)

colnames(fsim) <-entrez_ids$SYMBOL[c(1,3,4,5)]
rownames(fsim) <-entrez_ids$SYMBOL[c(1,3,4,5)]

for(i in 1:ncol(fsim))
  { fsim[i,i] <- NA } 
dat <- melt(fsim) 
dat <- dat[! is.na(dat$ value),] 
dat <- dat[,c( 1, 3)] 
head(dat)

dat.mean <- aggregate( value~Var1, dat, mean)
m<- dat.mean$value 
names(m) <- dat.mean $Var1 #按平均值给基因名排序
dat $Var1<- factor(dat $Var1, levels=names(sort(m)))
str(dat)

p <- ggplot(dat, aes(x = Var1, y = value, fill = factor(Var1))) +
  scale_fill_brewer(palette = "Set2", direction = 1) +  # 使用Set2调色板进行配色
  geom_boxplot() + coord_flip() +  # 坐标轴互换
  xlab("") + ylab("") + theme_bw()  # 设置坐标轴标签为空并应用白色主题

# 打印图形
print(p)
ggsave(filename = './Functional_correlation.pdf',plot  = p,width = 5,height = 5)

ggsave(filename = './Functional_correlation.png',plot  = p,width = 5,height =5)



func <-read.table('./GeneMANIA/genemania-functions.txt',fill = T,sep = '\t')
colnames(func) <-func[1,]
func <-func[-1,]
head(func)

func$log10_FDR <- -log10(as.numeric(func$FDR))

# Select the top 5 pathways based on log10_FDR
top5_func <- func[order(-func$log10_FDR), ][1:5, ]

# Plot
library(ggplot2)
ggplot(top5_func, aes(x = reorder(Function, log10_FDR), y = log10_FDR, fill = Function)) + 
  geom_bar(stat = 'identity') + 
  coord_flip() + 
  xlab('Function') + 
  ylab('-log10(FDR)') + 
  ggtitle('Top 5 Functions with -log10(FDR)') + 
  theme_minimal() + 
  scale_fill_brewer(palette = "Set3")+
  theme(legend.position = "none")
# Save as PDF
ggsave("03_top5_functions.pdf", plot = last_plot(), device = "pdf",width = 6,height = 3)

# Save as PNG
ggsave("03_top5_functions.png", plot = last_plot(), device = "png",width = 6,height = 3)

top5_func$Function
