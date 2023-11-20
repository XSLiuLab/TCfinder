
setwd("~/project/mcIdentify/data/")
GSE673_diff_path213 <- readRDS("~/project/mcIdentify/data/model_built/pathway_score_213/GSE673_diff_path213.rds")
library(data.table)
library(dplyr)


pathdata <- fread("./model_built/pathway_score_213/GSE673_pathway_score213.csv")
path_score1 <- pathdata

library(ggplot2)
library(ggpubr)
library(ggprism)
p1 <- ggplot(data=path_score1,aes(x=type,y=hsa04514,fill=factor(type)))+geom_boxplot(size=1,)+theme_prism()+
  stat_compare_means(label.x=1.2,size=5,method = "wilcox.test")+theme(legend.position = 'none')+
  labs(y="Pathway Score",title = "Cell adhesion molecules",x="")+
  theme(axis.title.x = element_blank())
p1

p2 <- ggplot(data=path_score1,aes(x=type,y=hsa04110,fill=factor(type)))+geom_boxplot(size=1,)+theme_prism()+
  stat_compare_means(label.x=1.2,size=5,method = "wilcox.test")+theme(legend.position = 'none')+
  labs(y="Pathway Score",title = "Cell cycle",x="")+
  theme(axis.title.x = element_blank())
p2

prow <- plot_grid(
  p1 ,
  p2,
  align = 'vh',
  labels = c(),
  hjust = -1,
  nrow = 1
)
prow


##Gene distribution


data1 <- fread("~/project/mcIdentify/data/model_built/GSE673_gene_distribution.csv")


library(ggpubr)
library(ggprism)
library(ggplot2)
ggplot(data=data1,aes(x=type,y=Freq))+
  geom_boxplot(size=1)+
  stat_compare_means(label.x=1.2,size=5,method = "wilcox.test")+
  theme_prism()+
  labs(y="Number of gene with expression")+
  scale_x_discrete(labels=c("Normal Cell","Tumor Cell"))+
  theme(axis.title.x = element_blank())


data1$type <- factor(data1$type,levels = c("normal","malignant"))
ggplot()+
  geom_density(data= GSE256_sample, alpha=0.8,adjust=1.5,aes(x=Freq,fill=type))+
  theme_prism()+
  labs(x="Number of gene with expression",y="Density")+
  scale_fill_manual(values = c("#F8766D","#00BFC4"))


ComplexHeatmap::Heatmap()


a <- as.data.frame(rnorm(20,mean = 2,sd = 1))
a$b <- c(1:20)
colnames(a) <- c("v1","v2")
ggplot(data = a,aes(x = v2,y=v1))+
  geom_point(size=3,color = "red")+
  labs(y="",x="")+
  theme_prism()



###heatmap
remove(list = ls())
setwd("~/project/mcIdentify/data/")
library(data.table)
library(dplyr)
data1 <- fread("./model_built/pathway_score_213/GSE673_pathway_score213.csv")
GSE673_diff_path213 <- readRDS("~/project/mcIdentify/data/model_built/pathway_score_213/GSE673_diff_path213.rds")

tumor <- data1 %>% filter(type=="malignant") %>% select(GSE673_diff_path213$hsa)
tumor1 <- tumor[1:500,]

normal <- data1 %>% filter(type=="normal") %>% select(GSE673_diff_path213$hsa)
normal1 <- normal[1:500,]

data2 <- rbind(tumor1,normal1)
library(scales)

data3 <- scale(data2)
data4 <- t(data3)


sample_group <- as.data.frame(c(rep("malignant",500),rep("normal",500)))
colnames(sample_group) <- "cluster"
library(ComplexHeatmap)
library(circlize)
col_fun = colorRamp2(c(-2, 0, 2), c("#00FF00", "#3B3B3B", "#EE0000"))
top_anno <- HeatmapAnnotation(Cluster = sample_group$cluster,
                              col = list(Cluster = c("malignant"= "#F8766D","normal"= "#00BFC4"),border = TRUE))
column_split = sample_group$cluster


library(ggprism)
ComplexHeatmap::Heatmap(data4,cluster_rows = T,cluster_columns = F,name = " ",
                        show_column_names = F,show_row_names = F,show_heatmap_legend = T,
                        col = col_fun,column_split = column_split,row_title = "Pathway")




remove(list = ls())
setwd("~/project/mcIdentify/data/")
library(data.table)
library(dplyr)
data1 <- fread("./model_built/pathway_score_213/model_result/GSE673_method_gene.csv")
data1$gene <- rep(1:6,4)
data1$method <- factor(data1$method,levels = c("mcIdentify","ikraus","SCINA","scMRMA"))
ggplot(data = data1, aes(x = gene, y = accuracy,  color = method, shape = method)) + 
  geom_point(size = 3) + 
  geom_smooth(size = 1.8) + 
  labs(x = " ", y = "Accuracy") + 
  ylim(0,1)+
  theme_prism()+
  scale_x_continuous(name = "Gene number", breaks = seq(1, 6, by = 1), 
                     labels = c("<500", "500~1000", "1000~1500", "1500~2000", "2000~2500",">2500"), limits = c(1, 6))




library(data.table)
library(dplyr)
data1 <- fread("./model_built/pathway_score_213/model_result/gene_sample_statistic.csv",header = T,data.table = F)

data1$type <- factor(data1$type,levels = c("normal","malignant"))
data1$gene <- rep(1:6,2)
library(ggplot2)

ggplot(data1, aes(x = gene, weight = value, fill = type))+
  geom_bar(position = "stack")+
  scale_fill_manual(values = c("#00BFC4","#F8766D"))+
  theme_prism()+
  scale_x_continuous(name = "Gene number", breaks = seq(1, 6, by = 1), 
                     labels = c("<500", "500~1000", "1000~1500", "1500~2000", "2000~2500",">2500"))+
  geom_text(aes(label = value1,y=value), 
            position = position_stack(vjust = 0.5), size = 5)









col_fun = colorRamp2(c(-2, 0, 2), c("red", "black", "blue"))
a <- matrix(data = rnorm(10000,mean = 0,sd = 1),nrow = 100,ncol = 100)

a[c(round(runif(50,min=1,max=100),0)),] = 0 

ComplexHeatmap::Heatmap(a,cluster_rows = F,cluster_columns = F,name = " ",
                        show_column_names = F,show_row_names = F,show_heatmap_legend = F,col = col_fun)
                        

?ComplexHeatmap::Heatmap()



remove(list = ls())
setwd("~/project/mcIdentify/data/")
library(data.table)
library(dplyr)
data1 <- fread("./model_built/pathway_score_213/GSE673_pathway_score213.csv")
GSE673_diff_path213 <- readRDS("~/project/mcIdentify/data/model_built/pathway_score_213/GSE673_diff_path213.rds")

tumor <- data1 %>% filter(type=="malignant") %>% select(GSE673_diff_path213$hsa[1:50])
tumor1 <- tumor[1:100,]

normal <- data1 %>% filter(type=="normal") %>% select(GSE673_diff_path213$hsa[1:50])
normal1 <- normal[1:100,]

data2 <- rbind(tumor1,normal1)
library(scales)

data3 <- scale(data2)
data4 <- t(data3)


sample_group <- as.data.frame(c(rep("malignant",100),rep("normal",100)))
colnames(sample_group) <- "cluster"
library(ComplexHeatmap)
library(circlize)
top_anno <- HeatmapAnnotation(Cluster = sample_group$cluster,
                              col = list(Cluster = c("malignant"= "#F8766D","normal"= "#00BFC4"),border = TRUE))


library(ggprism)
ComplexHeatmap::Heatmap(data4,cluster_rows = F,cluster_columns = F,name = " ",
                        show_column_names = F,show_row_names = F,show_heatmap_legend = F)


