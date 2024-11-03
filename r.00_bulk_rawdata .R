library(BiocManager)
# BiocManager::install("RTCGA")
# BiocManager::install("RTCGA.mRNA")
library(RTCGA)
library(RTCGAToolbox)
library(RTCGA.mRNA)
library(GEOquery)
library(dplyr)
library(tidyr)
library(R.utils)
library(data.table)
library(readxl)

#RTCGA包用于获取TCGA数据
setwd('/data/home/jiangminghui/project-0079/00_bulk_rawdata/')
library(readr)

con <- gzfile("./GSE233973_series_matrix.txt.gz", "rt")
# GSE70768 <- read.table(con,header = T,comment.char = "!",row.names = 1)
expr <- read.table(con, header = TRUE, comment.char = "!", row.names = 1, fill = TRUE)

# 关闭连接
close(con)

expr=as.matrix(expr)
samplenames <-colnames(expr)

expr = as.data.frame(expr) 
expr = data.frame(probe_id = rownames(expr),expr) 

anno <- read.table("./GPL21185-21174.txt", header = TRUE, sep = "\t", quote = "", fill = TRUE)

head(anno)

colnames(anno)

anno <-anno[-c(1,2),c('ID','GENE_SYMBOL')]

anno <-as.data.frame(anno)
source('./probe_annotation.R')
expr <- as.data.frame(expr)
expr[, -1] <- apply(expr[, -1], 2, as.numeric)
anno_unique <- anno %>% distinct(ID, .keep_all = TRUE)

symbol_exp =  probe_annotation(matrix = expr, annotate = anno_unique, mathod = 'mean')
rownames(symbol_exp) <-symbol_exp$symbol
symbol_exp <- symbol_exp[!is.na(rownames(symbol_exp)), ]

head(symbol_exp)
GSE233973_exp  <-symbol_exp
GSE233973_exp <- GSE233973_exp[,-1]

head(GSE233973_exp)
dim(GSE233973_exp)
colnames(GSE233973_exp)


GSM <- c("GSM7440033", "GSM7440034", "GSM7440035", "GSM7440036", "GSM7440037", "GSM7440038", "GSM7440039", "GSM7440040", "GSM7440041", "GSM7440042", "GSM7440043", "GSM7440044", "GSM7440045", "GSM7440046", "GSM7440047", "GSM7440048", "GSM7440049", "GSM7440050", "GSM7440051", "GSM7440052", "GSM7440053", "GSM7440054", "GSM7440055", "GSM7440056", "GSM7440057", "GSM7440058", "GSM7440059", "GSM7440060", "GSM7440061", "GSM7440062", "GSM7440063", "GSM7440064", "GSM7440065", "GSM7440066", "GSM7440067", "GSM7440068")
Sample <- c("JSA01", "NRA01", "NRA02", "SKA01", "SKA05", "SKA06", "TDA01", "TDA02", "TDA04", "TDA05", "TDA06", "TDA08", "TDA09", "TSA01", "SKH02", "SKH03", "SKH04", "SKH05", "TDH01", "TDH02", "TDH03", "TDH05", "TDH07", "TDH08", "TDH09", "TDH12", "TDH14", "JSN01", "SKN02", "SKN03", "TDN01", "TDN02", "TDN03", "TDN04", "TDN05", "TDN06")
sample_type <- data.frame(GSM, Sample)
# 找到不包含"TDH"的样本
use.col <- sample_type$Sample[!grepl("TDH|SKH", sample_type$Sample)]

# 获取这些样本在数据框中的位置
positions <- which(sample_type$Sample %in% use.col)
# 查看结果
print(positions)
GSE233973_exp <-GSE233973_exp[,sample_type$GSM[positions]]
# 计算每个基因的平均表达量
gene_means <- rowMeans(GSE233973_exp)

# 过滤掉平均表达量低于阈值的基因
GSE233973_exp <- GSE233973_exp[gene_means > 0, ]


saveRDS(GSE233973_exp ,'./GSE233973_exp.rds')

normal <-c('TDN','SKN','JSN')

sample_type <- sample_type %>%
  mutate(Group = ifelse(grepl("TDN|SKN|JSN", sample_type$Sample), 'normal', 'disease'))

saveRDS(sample_type,'./GSE233973_sample_type.rds')
head(colnames(GSE233973_exp))
group <- sample_type$Group[match(colnames(GSE233973_exp), sample_type$GSM)]
saveRDS(group ,'./GSE233973_group.rds')

file_path <- "./430408-revised-supplementary-tables.xls"
data <- read_excel(file_path, sheet = 1)
Macro_genes <-data$Symbol

saveRDS(Macro_genes,'./Macro_genes.rds')
