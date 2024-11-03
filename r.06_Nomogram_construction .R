rm(list = ls())
gc()
setwd('/data/home/jiangminghui/project-0079/06_Nomogram_construction/')

load('./Rdata_Nomogram_construction.RData')
library(dplyr)
library(tidyr)
library(ggpubr)
library(pROC)
library(rmda)
library(ggDCA)
library(rms)
expr <-readRDS('../00_bulk_rawdata/GSE233973_exp.rds')
expr<-log2(expr)
feature_gene <-readRDS('../05_machine_learning/03_feature_gene.rds')
group <-readRDS('../00_bulk_rawdata/GSE233973_group.rds')
data <- expr[feature_gene,] %>%t() %>%as.data.frame()
data$Group <- group
saveRDS(data,'./01_data.rds')

data <-readRDS('./01_data.rds')
write.csv(data ,'./01_data.csv',row.names = T)

dd <- datadist(data)
options(datadist = 'dd')
mylog<-lrm(Group~.,data=data,x=T,y=T,maxit=1000)

#画列线图
mynom<- nomogram(mylog, fun=plogis,fun.at=c(0.1,0.7,0.99),lp=F, funlabel="Risk of CAD")
cat(feature_gene)
# 构建逻辑回归模型
model <- lrm(Group ~ ADAM8 +FAM135A+ LCORL +PINK1 +SFI1, data = data)

# 绘制列线图
nomogram <- nomogram(model, fun = plogis, fun.at = c(0.1, 0.5, 0.9), lp = FALSE)
mynom<- nomogram(mylog, fun=plogis,fun.at=c(0.1,0.7,0.99),lp=F, funlabel="Risk of CAD")#"Risk of nonadherence"，可以改名，具体论文具体分析

plot(mynom,col.grid=c("tomato","darkcyan"))

pdf('./01_Nomogram_model.pdf',width = 10,height = 5)
plot(mynom,col.grid=c("tomato","darkcyan"))
dev.off()

png('./01_Nomogram_model.png',width = 10,height = 5,units = 'in',res = 300)
plot(mynom,col.grid=c("tomato","darkcyan"))
dev.off()



fit1 <- lrm(Group ~ ADAM8 +FAM135A+ LCORL +PINK1 +SFI1, data = data,x = T,y = T)
cali <-calibrate(fit1,B =1000)

pdf('./02_Calibration_Plot.pdf',width = 5,height = 5)

plot(cali,
     xlim = c(0,1),
     xlab = "Nomogram Predicted Probability",
     ylab = "Actual Probability",
     legend = FALSE,
     subtitles = FALSE)

abline(0,1,col = "#04198C",lty = 2,lwd = 2) 

lines(cali[,c("predy","calibrated.orig")], 
      type = "l",lwd = 2,col="#D91E1E",pch =16)

lines(cali[,c("predy","calibrated.corrected")],
      type = "l",lwd = 2,col="#4E7300",pch =16)

legend(0.55,0.35,  ##legend是绘制图例的函数，有兴趣的可以去深入了解这个函数
       c("Apparent","Ideal","Bias-corrected"), #表示曲线名称的集合
       lty = c(2,1,1), #lty表示线条类型，1代表实线，2至6都是虚线，虚的程度不一样
       lwd = c(2,1,1), #lwd表示宽度，以默认值的相对大小来表示宽度
       col = c("#04198C","#D91E1E","#4E7300"),  #给曲线添加颜色，对应上面c("Ax","Ix","Bx")
       bty = "n") # bty为o表示加边框，，注意是字母，不是数字0。bty可以取6种字符，分别为“o”、“l”、“7”、“c”、“u”、“]”
dev.off()

png('./02_Calibration_Plot.png',width = 5,height = 5,units = 'in',res = 300)

plot(cali,
     xlim = c(0,1),
     xlab = "Nomogram Predicted Probability",
     ylab = "Actual Probability",
     legend = FALSE,
     subtitles = FALSE)

abline(0,1,col = "#04198C",lty = 2,lwd = 2) 

lines(cali[,c("predy","calibrated.orig")], 
      type = "l",lwd = 2,col="#D91E1E",pch =16)

lines(cali[,c("predy","calibrated.corrected")],
      type = "l",lwd = 2,col="#4E7300",pch =16)

legend(0.55,0.35,  ##legend是绘制图例的函数，有兴趣的可以去深入了解这个函数
       c("Apparent","Ideal","Bias-corrected"), #表示曲线名称的集合
       lty = c(2,1,1), #lty表示线条类型，1代表实线，2至6都是虚线，虚的程度不一样
       lwd = c(2,1,1), #lwd表示宽度，以默认值的相对大小来表示宽度
       col = c("#04198C","#D91E1E","#4E7300"),  #给曲线添加颜色，对应上面c("Ax","Ix","Bx")
       bty = "n") # bty为o表示加边框，，注意是字母，不是数字0。bty可以取6种字符，分别为“o”、“l”、“7”、“c”、“u”、“]”
dev.off()


form.bestglm<-as.formula(Group ~ ADAM8 +FAM135A+ LCORL +PINK1 +SFI1)

data$Group <-ifelse(data$Group  =='disease',1,0)

DCA.1<- decision_curve(formula=form.bestglm,
                       family = binomial(link ='logit'),
                       thresholds= seq(0,1, by = 0.01),
                       confidence.intervals =0.95,
                       study.design = 'cohort',
                       data = data)

DCA.1$derived.data
head(DCA.1$derived.data[,c("thresholds","NB","sNB","cost.benefit.ratio")])

plot_decision_curve(DCA.1,
                    curve.names= c("Model Bestglm"),
                    xlim=c(0,0.8),
                    cost.benefit.axis =TRUE,
                    col = "#E64B35B2",
                    confidence.intervals =FALSE,
                    standardize = FALSE)

png("03_decision_curve.png", width = 6, height = 6,res = 300,units = 'in')
plot_decision_curve(DCA.1,
                    curve.names= c("Model Bestglm"),
                    xlim=c(0,0.8),
                    cost.benefit.axis =TRUE,
                    col = "#E64B35B2",
                    confidence.intervals =FALSE,
                    standardize = FALSE)
dev.off()

# Save as PDF
pdf("03_decision_curve.pdf", width = 6, height = 6)
plot_decision_curve(DCA.1,
                    curve.names= c("Model Bestglm"),
                    xlim=c(0,0.8),
                    cost.benefit.axis =TRUE,
                    col = "#E64B35B2",
                    confidence.intervals =FALSE,
                    standardize = FALSE)
dev.off()


set.seed(123)
fit2 <- decision_curve(Group ~ ADAM8 +FAM135A+ LCORL +PINK1 +SFI1,
                       data = data, 
                       bootstraps = 50
)

plot_decision_curve(list(DCA.1, fit2),
                    curve.names = c("fit1", "fit2"), 
                    xlim = c(0, 1), # 可以设置x轴范围
                    legend.position = "topright", # 图例位置,
                    col = c("#91A646","#D97904"), # 自定义颜色
                    confidence.intervals = "none",
                    lty = c(1,2), # 线型，注意顺序
                    lwd = c(3,2,2,1) # 注意顺序，先是自己的2个模型，然后是All,然后是None
)

# Create the plot
 

png("04_clinical_impact.png", width = 5, height = 5,units = 'in',res = 300)
plot_clinical_impact(
  DCA.1,
  population.size = 1000,
  cost.benefit.axis = TRUE,
  n.cost.benefits = 8,
  col = c('#91A646','#D97904'),
  confidence.intervals = TRUE,
  ylim = c(0,1000),
  legend.position = "topright"
)
dev.off()

# Save as PDF
pdf("04_clinical_impact.pdf", width = 5, height = 5)
plot_clinical_impact(
  DCA.1,
  population.size = 1000,
  cost.benefit.axis = TRUE,
  n.cost.benefits = 8,
  col = c('#91A646','#D97904'),
  confidence.intervals = TRUE,
  ylim = c(0,1000),
  legend.position = "topright"
)
dev.off()
# 
# data1 <-list()
# data1$Group<- group
# 
# 
# plot_combined_roc <- function(data, biomarkers, exp_data) {
#   plot(NULL, xlim = c(0, 1), ylim = c(0, 1), xlab = "1 - Specificity", ylab = "Sensitivity", main = " ")
#   colors <- rainbow(length(biomarkers))
#   
#   for (i in seq_along(biomarkers)) {
#     biomarker <- biomarkers[i]
#     data$Biomarker_exp <- exp_data[biomarker, ] %>% as.numeric()
#     
#     roc_result <- roc(data$Group, data$Biomarker_exp)
#     auc_value <- auc(roc_result)
#     
#     plot(roc_result, col = colors[i], add = TRUE)
#     text(0.8, 0.4 - 0.05 * i, paste(biomarker, "AUC =", round(auc_value, 2)), col = colors[i])
#   }
#   
#   legend("bottomleft", legend = biomarkers, col = colors, lwd = 2)
# }
# 
# pdf('./04_AUC_curve.pdf',width = 5,height = 5)
# plot_combined_roc(data1, feature_gene, expr)
# dev.off()
# 
# png('./04_AUC_curve.png',width = 5,height = 5,units = 'in',res = 300)
# plot_combined_roc(data1, feature_gene, expr)
# dev.off()


# 获取当前目录的路径
current_dir <- getwd()

# 获取当前目录的文件夹名字
dir_name <- basename(current_dir)

# 保存当前环境到文件
save.image(file = paste0(current_dir, "/", "Rdata_",dir_name, ".RData"))





