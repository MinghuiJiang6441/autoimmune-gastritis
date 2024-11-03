rm(list = ls())
gc()
library(readxl)
library(readr)
library(jsonlite)
library(dplyr)
library(stringr)
library(ggalluvial)
setwd('/data/home/jiangminghui/project-0079/11_TF_gene_miRNA/')
feature_gene <-readRDS('../05_machine_learning/03_feature_gene.rds')

# 获取当前工作目录
current_dir <- getwd()

# 列出当前目录下的所有文件和文件夹
file_list <- list.files(current_dir, recursive = TRUE)

# 遍历所有文件并执行操作
for (file in file_list) {
  file_path <- file.path(current_dir, file)

    print(file_path)
}

ADAM8_database1 <-read_excel( './ADAM8_8.0_ENST00000415217.3_predicted_targeting_details.xlsx',skip = 3) %>%as.data.frame()
colnames(ADAM8_database1)[1] <- 'mirna'
ADAM8_database1 <-ADAM8_database1$mirna %>% as.character()

file_list <- list.files(current_dir, pattern = "*.xlsx$", recursive = TRUE, full.names = TRUE)

# 定义函数处理每个文件
process_file <- function(file_path) {
  data <- read_excel(file_path, skip = 3) %>% as.data.frame()
  colnames(data)[1] <- 'mirna'
  mirna <- data$mirna %>% as.character()
  return(mirna)
}

# 遍历所有文件并处理
file_contents <- lapply(file_list, process_file)
names(file_contents)  <-c('ADAM8','FAM135A','LCORL','PINK1','SFI1')


gene2mir_data <- read_csv('./gene2mir.csv')
head(gene2mir_data)
gene2mir_data <-split(gene2mir_data$ID,gene2mir_data$Target)

database3 <-list()
data <- read.table('./ENCORI_hg38_CLIP-seq_miRNA-target_all_FAM135A.txt')
data <-data[-1,]
database3[[1]] <-data

data <- read.table('./ENCORI_hg38_CLIP-seq_miRNA-target_all_LCORL.txt')
data <-data[-1,]
head(data)
database3[[2]] <-data

data <- read.table('./ENCORI_hg38_CLIP-seq_miRNA-target_all_PINK1.txt')
data <-data[-1,]
head(data)
database3[[3]] <-data

names(database3) <-c('FAM135A','LCORL','PINK1')
combined_data <- do.call(rbind, database3)
database3_l <-split(combined_data$V2,combined_data$V4)



names(file_contents)
names(gene2mir_data)
names(database3_l)


genes <- union(names(file_contents), union(names(gene2mir_data), names(database3_l)))

# 创建一个函数来计算交集
calculate_intersect <- function(gene, list1, list2, list3) {
  sets <- list()
  
  if (gene %in% names(list1)) {
    sets[[1]] <- list1[[gene]]
  } else {
    sets[[1]] <- NULL
  }
  
  if (gene %in% names(list2)) {
    sets[[2]] <- list2[[gene]]
  } else {
    sets[[2]] <- NULL
  }
  
  if (gene %in% names(list3)) {
    sets[[3]] <- list3[[gene]]
  } else {
    sets[[3]] <- NULL
  }
  
  # 移除 NULL 集合
  sets <- Filter(Negate(is.null), sets)
  
  if (length(sets) < 2) {
    message("Gene ", gene, " has less than two sets available for intersection.")
    return(NULL)
  }
  
  # 计算交集
  intersect_mirna <- Reduce(intersect, sets)
  return(intersect_mirna)
}

# 遍历所有基因并计算交集
intersect_results <- list()
for (gene in genes) {
  intersect_results[[gene]] <- calculate_intersect(gene, file_contents, gene2mir_data, database3_l)
}

# 打印交集结果
print(intersect_results)
saveRDS(intersect_results ,'./01_intersect_results.rds')


cat(feature_gene ,sep = '\n')

# csv_data <- read_csv('./node.csv')
# head(csv_data)
json_data <- fromJSON('./LCORL_PINK1_TF.json')
head(json_data)

edge_names <- json_data$elements$edges$data$name

# 去掉 (pp) 并拆分成两列
edge_df1 <- data.frame(
  source = str_replace(edge_names, " \\(pp\\) .*", ""),
  target = str_replace(edge_names, ".* \\(pp\\) ", "")
)
head(edge_df1)

json_data <- fromJSON('./SFI1_ADAM.json')
head(json_data)

edge_names <- json_data$elements$edges$data$name

# 去掉 (pp) 并拆分成两列
edge_df2 <- data.frame(
  source = str_replace(edge_names, " \\(pp\\) .*", ""),
  target = str_replace(edge_names, ".* \\(pp\\) ", "")
)
head(edge_df2)

json_data <- fromJSON('./FAM135A.json')
head(json_data)

edge_names <- json_data$elements$edges$data$name

# 去掉 (pp) 并拆分成两列
edge_df3 <- data.frame(
  source = str_replace(edge_names, " \\(pp\\) .*", ""),
  target = str_replace(edge_names, ".* \\(pp\\) ", "")
)
head(edge_df3)
edge_df <-rbind(edge_df1,edge_df2,edge_df3)

library(dplyr)

# 将intersect_results转化为数据框格式
mirna_intersections <- do.call(rbind, lapply(names(intersect_results), function(gene) {
  if (length(intersect_results[[gene]]) > 0) {
    data.frame(Gene = gene, miRNA = intersect_results[[gene]])
  } else {
    NULL
  }
}))



combined_df <- merge(mirna_intersections, edge_df, by.x = "Gene", by.y = "source", all = TRUE)
saveRDS(combined_df,'./combined_df.rds')
# mytheme1<-theme_bw() +
#   theme(panel.grid =element_blank()) +
#   theme(panel.border = element_blank()) +
#   theme(axis.line = element_blank(),axis.ticks =element_blank(),axis.text.y =element_blank())
# 
# gene_colors <- c(
#   "ADAM8" = "#FF5733", 
#   "FAM135A" = "#33FF57", 
#   "LCORL" = "#3357FF", 
#   "PINK1" = "#F333FF", 
#   "SFI1" = "#FFD433"
# )
# ggplot(combined_df, aes(axis1 = target, axis2 = Gene, axis3 = miRNA, y = 1,fill = Gene)) +
#   geom_alluvium(aes(fill = Gene)) +
#   geom_stratum() +
#   geom_text(stat = "stratum", aes(label = after_stat(stratum)), size = 3, color = "black") +
#   scale_x_discrete(limits = c("Target", "Gene", "miRNA"), expand = c(.05, .05)) +
#   scale_fill_manual(values = gene_colors) +
#   theme_void() +
#   theme(
#     axis.text.y = element_blank(),
#     legend.title = element_text(size = 10),
#     legend.text = element_text(size = 8),
#     plot.title = element_text(size = 6, face = "bold")
#   ) +
#   ggtitle("Sankey Diagram for Gene-MiRNA Interactions")

# 合并数据框
combined_df <- merge(mirna_intersections, edge_df, by.x = "Gene", by.y = "source", all = TRUE)

# 计算每个 target 和 miRNA 的出现次数
target_count <- table(combined_df$target)
miRNA_count <- table(combined_df$miRNA)

# 保留出现次数大于1的 target 和 miRNA
 filtered_df <- combined_df[combined_df$target %in% names(target_count[target_count > 1]) & 
                              combined_df$miRNA %in% names(miRNA_count[miRNA_count > 1]), ]

filtered_df$Freq <-rep(1)
df_lodes <- to_lodes_form(
  filtered_df,
  key = "x",
  value = "stratum",
  id = "alluvium",
  axes = c(2, 1, 3)  # 调整列顺序，确保target在第一列，Gene在第二列，miRNA在第三列
)
mycol3=colorRampPalette(c("#00abef","#64b036","#ffe743","#64b036","#00abef"))(117)

p <-ggplot(df_lodes,aes(x = x, stratum =stratum, alluvium = alluvium,
                    fill = stratum, label = stratum)) +
  scale_x_discrete(expand = c(0, 0)) +
  geom_flow(width = 0.2, knot.pos = 0.1) +
  geom_stratum(alpha = .9,color="grey20",width = 1/7) +
  geom_text(stat = "stratum", size =1.5,color="black") +
  scale_fill_manual(values = mycol3) +
  xlab("") + ylab("") +
  theme_bw() +
  theme(panel.grid =element_blank()) +
  theme(panel.border = element_blank()) +
  theme(axis.line = element_blank(),axis.ticks =element_blank(),axis.text.y =element_blank())+
  guides(fill = FALSE)
p
# 保存为PDF格式
ggsave(filename = "sankey_diagram.pdf", plot = p, device = "pdf", width = 4, height = 8)

# 保存为PNG格式
ggsave(filename = "sankey_diagram.png", plot = p, device = "png", width = 4, height = 8, dpi = 300)
