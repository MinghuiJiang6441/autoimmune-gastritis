rm(list = ls())
gc()

setwd('/data/home/jiangminghui/project-0079/05_machine_learning/')
load('./rdata_machine_learning.RData')
library(e1071)
library(caret)
library(ggplot2)
library(cowplot)
library(glmnet)
library(ggvenn)

expr <-readRDS('../00_bulk_rawdata/GSE233973_exp.rds')
expr <- log2(expr)

Overlapgenes <- readRDS('../01_DEGs/02_Overlapgenes.rds')
# res <-readRDS('../01_DEGs/01_DEGs_res.rds')
# 添加显著性标记
# res$significant <- ifelse(res$adj.P.Val < 0.05 & res$logFC > 0.5, "Upregulated",
#                           ifelse(res$adj.P.Val < 0.05 & res$logFC < -0.5, "Downregulated", "Not Significant"))
# 
# a <-readRDS('../00_bulk_rawdata/Macro_genes.rds')
# 
# b <- subset(res,significant !="Not Significant") %>%rownames()
# 
# 
# Overlapgenes <-intersect(a,b)
group <-readRDS('../00_bulk_rawdata/GSE233973_group.rds')
gene_data <- data.frame(expr[Overlapgenes,] %>% t())
gene_data$labels <-group
gene_data$labels <- factor(gene_data$labels, levels = c('disease', 'normal'))

set.seed(21) # 设置种子
control <- rfeControl(functions = caretFuncs, method = "cv", number = 5)

num <- ncol(gene_data)-1
results <- rfe(x = gene_data[, 1:num], # 除去最后一列，其余列均为预测变量（也就是hubgene的表达量）
               y = gene_data$labels, # 分组信息
               sizes = c(1:num),
               rfeControl = control,
               method = "svmRadial"
)
saveRDS(results,'./01_svmProfile.rds')
results<-readRDS('./01_svmProfile.rds')
svmrfe_result <- data.frame(symbol = predictors(results))

p1 <- plot(results, type = c("o"), xgap.axis = 1,lwd =6)

# 使用plot_grid函数组合图形
combined_plot <- plot_grid(p1, labels = "AUTO")

# 添加标题和主题元素
combined_plot <- combined_plot +
  draw_label("SVM_RFE_analyse", x = 0.5, y = 1, hjust = 0.5, vjust = 1, size = 25, fontface = "bold", color = "black") +
  theme(plot.title = element_text(hjust = 0.5, colour = "black", face = "bold", size = 25),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.text = element_blank(),
        legend.title = element_blank(),
        legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  annotate("text", x = 8, y = 0.92, label = "n=8, Accuracy = 0.92", size = 5, hjust = 0.75, vjust = 0.75)

# 显示组合图形
print(combined_plot)

ggsave('./01_SVM_RFE_analyse.pdf',combined_plot,width = 7,height = 7)
ggsave('./01_SVM_RFE_analyse.png',combined_plot,width = 7,height = 7)


x <- as.matrix(gene_data[,-ncol(gene_data)])
y <- ifelse(group =='disease',1,0)

set.seed(111)#不要改！！！！！！！

fit <-glmnet(x,y,alpha = 1)

pdf('./02_Lassoplot_1.pdf',width = 6,height = 5)
plot(fit,xvar="lambda",label=T)
dev.off()

png('./02_Lassoplot_1.png',width = 6,height = 5,units = 'in',res = 300)
plot(fit,xvar="lambda",label=T)
dev.off()

# 执行LASSO回归分析
set.seed(111)#不要改！！！！！！！
cv_fit <- cv.glmnet(x, y, alpha = 1, nfolds = 10)



# 获取与最小偏似然偏差相对应的惩罚参数值（λ）
best_lambda <- cv_fit$lambda.min
set.seed(11333)#不要改！！！！！！！
# 使用最佳λ值重新拟合模型
lasso_model <- glmnet(x, y, alpha = 1, lambda = best_lambda)
pdf('./02_Lassoplot.pdf',width = 6,height = 5)
plot(cv_fit, xvar = "lambda", label = TRUE)
abline(v = log(best_lambda), col = " red", lty = 2) 
dev.off()

png('./02_Lassoplot.png',width = 6,height = 5,units = 'in',res = 300)
plot(cv_fit, xvar = "lambda", label = TRUE)
abline(v = log(best_lambda), col = " red", lty = 2) 
dev.off()


# 获取β系数矩阵
coef_matrix <- as.matrix(coef(lasso_model))

# 获取β系数不为0的特征基因
selected_genes <- rownames(coef_matrix)[coef_matrix != 0]
selected_genes <- selected_genes[selected_genes != "(Intercept)"]

library(Boruta)
# set.seed(123) # 设置随机种子以确保结果可重复
# boruta_result <- Boruta(labels ~ ., data = gene_data, doTrace = 2, pValue = 0.01)
# 
# print(boruta_result)
# 
# final_features <- getSelectedAttributes(boruta_result, withTentative = TRUE)
# print(final_features)
# 
# data_selected <- gene_data[, c(final_features, "labels")]
# 
# dim(data_selected)

# library(randomForest)
# 
# # 构建随机森林模型
# set.seed(123)
# model <- randomForest(labels ~ ., data = gene_data, importance = TRUE)
# 
# # 打印模型结果
# print(model)
# 
# # 查看特征重要性
# importance(model)
# 
# rf <- randomForest(gene_data, gene_data$labels)
# followerDF <- data.frame(Real_Follower=gene_data$labels, Predicted_Follower=predict(rf, newdata=gene_data))
# library(ggplot2)
# sp_scatterplot(followerDF, xvariable = "Real_Follower", yvariable = "Predicted_Follower",
#                smooth_method = "auto") + coord_fixed(1)



set.seed(1)

boruta <- Boruta(x=gene_data[,-ncol(gene_data)], y=gene_data$labels, pValue=0.01, mcAdj=T, 
                 maxRuns=300)

boruta
table(boruta$finalDecision)

boruta.imp <- function(x){
  imp <- reshape2::melt(x$ImpHistory, na.rm=T)[,-1]
  colnames(imp) <- c("Variable","Importance")
  imp <- imp[is.finite(imp$Importance),]
  
  variableGrp <- data.frame(Variable=names(x$finalDecision), 
                            finalDecision=x$finalDecision)
  
  showGrp <- data.frame(Variable=c("shadowMax", "shadowMean", "shadowMin"),
                        finalDecision=c("shadowMax", "shadowMean", "shadowMin"))
  
  variableGrp <- rbind(variableGrp, showGrp)
  
  boruta.variable.imp <- merge(imp, variableGrp, all.x=T)
  
  sortedVariable <- boruta.variable.imp %>% group_by(Variable) %>% 
    summarise(median=median(Importance)) %>% arrange(median)
  sortedVariable <- as.vector(sortedVariable$Variable)
  
  
  boruta.variable.imp$Variable <- factor(boruta.variable.imp$Variable, levels=sortedVariable)
  
  invisible(boruta.variable.imp)
}

boruta.variable.imp <- boruta.imp(boruta)

gene_con <- subset(boruta.variable.imp,boruta.variable.imp$finalDecision == 'Confirmed')
boruta_signif <-unique(gene_con$Variable) %>% as.character()


head(boruta.variable.imp)

library(ImageGP)

plot <- sp_boxplot(boruta.variable.imp, melted = TRUE, xvariable = "Variable", yvariable = "Importance",
                   legend_variable = "finalDecision", legend_variable_order = c("shadowMax", "shadowMean", "shadowMin", "Confirmed"),
                   xtics_angle = 90, coordinate_flip = TRUE)

# Save the plot as PNG
ggsave("03_sp_boxplot.png", plot, width = 5, height = 5, dpi = 300)

# Save the plot as PDF
ggsave("03_sp_boxplot.pdf", plot, width = 5, height = 5)
# 
# 
# boruta.finalVarsWithTentative <- data.frame(Item=getSelectedAttributes(boruta, withTentative = T), Type="Boruta_with_tentative")
# variableFactor <- rev(levels(boruta.variable.imp$Variable))
# labels <- as.vector(gene_data$labels)
# sp_scatterplot(gene_data, xvariable =  labels, yvariable = variableFactor[1], smooth_method = "auto")

# 打印选中的特征基因
print(selected_genes)

feature_gene <-intersect(boruta_signif,selected_genes)
feature_gene

saveRDS(feature_gene,'./03_feature_gene.rds')


write.csv(feature_gene,file = './03_feature_gene.csv')

feature_gene <-readRDS('./03_feature_gene.rds')

venn <- list('boruta_signif' =boruta_signif,
             'Lasso' = selected_genes
)

# Draw the Venn diagram
venn_plot <- ggvenn(
  venn,
  c("boruta_signif", "Lasso"),
  fill_color = c("#E69F00", "#56B4E9"),
  stroke_color = "black",
  stroke_size = 1,
  set_name_color = c("#E69F00", "#56B4E9"),
  text_color = "black",
  text_size = 6
)

# Save the Venn diagram as PNG
ggsave('./04_Venn_Plot.png', venn_plot, width = 6, height = 6, dpi = 300)

# Save the Venn diagram as PDF
ggsave('./04_Venn_Plot.pdf', venn_plot, width = 6, height = 6)



# 获取当前目录的路径
current_dir <- getwd()

# 获取当前目录的文件夹名字
dir_name <- basename(current_dir)

# 保存当前环境到文件
save.image(file = paste0(current_dir, "/", 'rdata_',dir_name, ".RData"))

