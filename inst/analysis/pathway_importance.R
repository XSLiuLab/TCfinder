
setwd("~/project/mcIdentify/data/")

library(data.table)
library(dplyr)



path_impor1 <- fread("./model_built/pathway_score_213/pathway_importance/path_impor_GSE673_model15.csv")
path_score1 <- fread("./model_built/pathway_score_213/GSE673_pathway_score213.csv")
path_impor1$hsa <- colnames(path_score1)

path_impor1 <- path_impor1[-1,-1]
path_impor1$number <- c(1:213)

# top10 <- path_impor1 %>% arrange(-V2) %>% filter(V2 > 0.059)

plot(path_impor1$V2,ylim = c(0.945,0.985),col = "blue", pch = 19, cex = 1)


path_impor1 <- path_impor1 %>% dplyr::mutate(fac = case_when(V2 < 0.979 ~ "A", TRUE ~ "B"))
path_impor1$fac <- as.factor(path_impor1$fac)

path_impor1 <- left_join(path_impor1,GSE673_diff_path213)


library(ggplot2)
library(ggpubr)
library(ggprism)
library(ggrepel)
ggplot(path_impor1, aes(x=number, y=V2, color=factor(fac))) + 
  geom_point(size = 3,)+
  theme_prism(border = T)+
  labs(y="Accuracy of the model", x  = "Pathway")+
  ylim(0.955,0.984)+
  xlim(0,214)+
  scale_color_manual(values = c("#DC0000FF",'#0072B5FF'))+
  theme(legend.position = 'none')+
  geom_text_repel(
    data = subset(path_impor1, path_impor1$V2 < 0.979),
    aes(label = pathway_id),
    size = 5,
    box.padding = unit(1, "lines"),
    point.padding = unit(1, "lines"), segment.color = "black", show.legend = FALSE )+
  geom_hline(aes(yintercept=0.98),linetype=5,col="black")
  
  


library(data.table)
library(dplyr)

path_impor2 <- fread("./model_built/pathway_score_213/pathway_importance/loss_result_GSE673_model15.csv")
path_score2 <- fread("./model_built/pathway_score_213/GSE673_pathway_score213.csv")
path_impor2$hsa <- colnames(path_score2)

path_impor2 <- path_impor2[-1,-1]
path_impor2$number <- c(1:213)

path_impor2 <- left_join(GSE673_diff_path213,path_impor2)

draw_data <- path_impor2 %>% filter(hsa %in% top10$hsa) %>% select(-pvalue) %>% select(-number)


data1 <- melt(draw_data)
data1$hsa <- factor(data1$hsa,levels = top10$hsa)


library(ggpubr)
library(ggprism)
library(ggplot2)
library(cowplot)
ggplot(data=data1,aes(x=hsa,y=value,fill = pathway_id))+
  geom_boxplot(size=1, draw_quantiles = c(0.5))+
  theme_prism(border = T)+theme(legend.position = 'none')+
  labs(y="Loss of the model",title = " ")+
  theme(axis.title.x = element_blank())+
  # geom_hline(aes(yintercept=0.9),linetype=5,col="red")+
  # geom_hline(aes(yintercept=0.8),linetype=5,col="red")+
  # scale_y_continuous(breaks=c(0,0.5,0.6,0.7,0.8,0.9,1))+
  scale_fill_manual(values = c("#F27970", "#BB9727","#54B345","#32B897",
                               "#05B9E2", "#8983BF","#C76DA2","#F27970",
                               "#BB9727","#54B345"))





library(ggpubr)
library(ggprism)
library(ggplot2)
library(cowplot)
ggplot(data=path_score1,aes(x=type,y=hsa05416,fill=factor(type)))+
  geom_boxplot(size=1,)+
  stat_compare_means(label.x=1.2,size=5,method = "wilcox.test")+
  theme_prism()+theme(legend.position = 'none')+
  labs(y="Pathway Score",title = "Viral myocarditis")+
  theme(axis.title.x = element_blank())




p1 <- ggplot(data=path_score1,aes(x=type,y=hsa00190,fill=factor(type)))+geom_boxplot(size=1,)+theme_prism()+
  stat_compare_means(label.x=1.2,size=5,method = "wilcox.test")+theme(legend.position = 'none')+
  labs(y="Pathway Score",title = "Oxidative phosphorylation")+
  theme(axis.title.x = element_blank())

p2 <- ggplot(data=path_score1,aes(x=type,y=hsa04612,fill=factor(type)))+geom_boxplot(size=1,)+theme_prism()+
  stat_compare_means(label.x=1.2,size=5,method = "wilcox.test")+theme(legend.position = 'none')+
  labs(y="Pathway Score",title = "Antigen processing and presentation")+
  theme(axis.title.x = element_blank())

p3 <- ggplot(data=path_score1,aes(x=type,y=hsa04940,fill=factor(type)))+geom_boxplot(size=1,)+theme_prism()+
  stat_compare_means(label.x=1.2,size=5,method = "wilcox.test")+theme(legend.position = 'none')+
  labs(y="Pathway Score",title = "Type I diabetes mellitus")+
  theme(axis.title.x = element_blank())

p4 <- ggplot(data=path_score1,aes(x=type,y=hsa05416,fill=factor(type)))+geom_boxplot(size=1,)+theme_prism()+
  stat_compare_means(label.x=1.2,size=5,method = "wilcox.test")+theme(legend.position = 'none')+
  labs(y="Pathway Score",title = "Viral myocarditis")+
  theme(axis.title.x = element_blank())



ggdraw() +     
  draw_plot(p3, 0, 0, 0.5, 0.5) +  
  draw_plot(p4, 0.5, 0, 0.5, 0.5) +  
  draw_plot(p1, 0, 0.5, 0.5, 0.5) +
  draw_plot(p2, 0.5, 0.5, 0.5, 0.5)



library(ggpubr)
library(ggprism)
library(ggplot2)
ggplot(data=path_score1,aes(x=type,y=hsa04940))+
  geom_boxplot(size=1,)+
  stat_compare_means(label.x=1.2,size=5,method = "wilcox.test")+
  theme_prism()+
  labs(y="Pathway Score")+
  theme(axis.title.x = element_blank())




######pathway importance 
path_impor <- fread("./model_built/pathway_score_213/pathway_importance/path_impor_GSE530_model15.csv")
path_score <- fread("./model_built/pathway_score_213/GSE530_pathway_score213.csv")
GSE673_diff_path213 <- readRDS("~/project/mcIdentify/data/model_built/pathway_score_213/GSE673_diff_path213.rds")
path_impor$hsa <- colnames(path_score)

path_impor <- path_impor[-1,-1]
path_impor$number <- c(1:213)

# plot(path_impor$V2,ylim = c(0.105,0.175),col = "blue", pch = 19, cex = 1)

path_impor <- path_impor %>% dplyr::mutate(fac = case_when(V2 > 0.11 ~ "A", TRUE ~ "B"))
path_impor$fac <- as.factor(path_impor$fac)

path_impor <- left_join(path_impor,GSE673_diff_path213)

library(ggplot2)
library(ggpubr)
library(ggprism)
library(ggrepel)
ggplot(path_impor, aes(x=number, y=V2, color=factor(fac))) + 
  geom_point(size = 3,)+
  theme_prism(border = T)+
  labs(y="Loss of the model", x  = "Pathway")+
  ylim(0.089,0.136)+
  scale_color_manual(values = c('red','blue'))+
  theme(legend.position = 'none')+
  geom_text_repel(
  data = subset(path_impor, path_impor$V2 > 0.11),
  aes(label = pathway_id),
  size = 4,
  box.padding = unit(1.2, "lines"),
  point.padding = unit(1, "lines"), segment.color = "black", show.legend = FALSE )





library(ggpubr)
library(ggprism)
library(ggplot2)
library(cowplot)
ggplot(data=path_score,aes(x=type,y=hsa04380,fill=factor(type)))+
  geom_boxplot(size=1,)+
  stat_compare_means(label.x=1.2,size=5,method = "wilcox.test")+
  theme_prism()+theme(legend.position = 'none')+
  labs(y="Pathway Score",title = "Osteoclast differentiation")+
  theme(axis.title.x = element_blank())

                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                          


p1 <- ggplot(data=path_score,aes(x=type,y=hsa04380,fill=factor(type)))+geom_boxplot(size=1,)+theme_prism()+
  stat_compare_means(label.x=1.2,size=5,method = "wilcox.test")+theme(legend.position = 'none')+
  labs(y="Pathway Score",title = "Osteoclast differentiation")+
  theme(axis.title.x = element_blank())

p2 <- ggplot(data=path_score,aes(x=type,y= hsa04940,fill=factor(type)))+geom_boxplot(size=1,)+theme_prism()+
  stat_compare_means(label.x=1.2,size=5,method = "wilcox.test")+theme(legend.position = 'none')+
  labs(y="Pathway Score",title = "Type I diabetes mellitus")+
  theme(axis.title.x = element_blank())

p3 <- ggplot(data=path_score,aes(x=type,y=hsa04650,fill=factor(type)))+geom_boxplot(size=1,)+theme_prism()+
  stat_compare_means(label.x=1.2,size=5,method = "wilcox.test")+theme(legend.position = 'none')+
  labs(y="Pathway Score",title = "Natural killer cell mediated cytotoxicity")+
  theme(axis.title.x = element_blank())

p4 <- ggplot(data=path_score,aes(x=type,y=hsa04978,fill=factor(type)))+geom_boxplot(size=1,)+theme_prism()+
  stat_compare_means(label.x=1.2,size=5,method = "wilcox.test")+theme(legend.position = 'none')+
  labs(y="Pathway Score",title = "Mineral absorption")+
  theme(axis.title.x = element_blank())

p5 <- ggplot(data=path_score,aes(x=type,y=hsa05322,fill=factor(type)))+geom_boxplot(size=1,)+theme_prism()+
  stat_compare_means(label.x=1.2,size=5,method = "wilcox.test")+theme(legend.position = 'none')+
  labs(y="Pathway Score",title = "Systemic lupus erythematosus")+
  theme(axis.title.x = element_blank())

p6 <- ggplot(data=path_score,aes(x=type,y=hsa05208,fill=factor(type)))+geom_boxplot(size=1,)+theme_prism()+
  stat_compare_means(label.x=1.2,size=5,method = "wilcox.test")+theme(legend.position = 'none')+
  labs(y="Pathway Score",title = "Chemical carcinogenesis - reactive oxygen species")+
  theme(axis.title.x = element_blank())


ggdraw() +     
  draw_plot(p3, 0, 0, 0.33, 0.5) +  
  draw_plot(p2, 0.33, 0, 0.33, 0.5) +  
  draw_plot(p1, 0.66, 0, 0.33, 0.5) +
  draw_plot(p6, 0, 0.5, 0.33, 0.5) +
  draw_plot(p4, 0.33, 0.5, 0.33, 0.5) +
  draw_plot(p5, 0.66, 0.5, 0.33, 0.5)




gene_ <- KEGG_pathway_gene %>% filter(hsa %in% path_impor1[path_impor1$V2>0.058,]$hsa)
sort(table(gene_$gene_id))



gene_list <- KEGG_pathway_gene %>% filter(hsa %in% path_impor[path_impor$V2 > 0.135,]$hsa)






###heatmap pathy
pathway_data <- fread("./model_built/pathway_score_213/GSE673_pathway_score213.csv")
GSE673_diff_path213 <- readRDS("~/project/mcIdentify/data/model_built/pathway_score_213/GSE673_diff_path213.rds")

heatpathay <- GSE673_diff_path213 %>% filter(hsa %in% c("hsa00190","hsa04612","hsa04940","hsa05416"))

heatmap_data <- pathway_data %>% select(type,heatpathay$hsa) %>% arrange(type)
heatmap_data1 <-  heatmap_data %>% select(-type)

heatmap_data2 <- scale(heatmap_data1)
heatmap_data2 <- t(heatmap_data2)

heatmap_data3 <- heatmap_data2[,c(1:2000,33001:35000)]



tumor_sample <- heatmap_data %>% filter(type =="malignant") %>% arrange(-hsa00190)
tumor_sample1 <- tumor_sample[1:2000,]

normal_sample <- heatmap_data %>% filter(type =="normal") %>% arrange(-hsa04612)
normal_sample1 <- normal_sample[1:2000,]

data1 <- rbind(tumor_sample1,normal_sample1)
data2 <- data1[,-1]
data3 <- scale(data2)
# data3 <- data2
data4 <- t(data3)

library(ComplexHeatmap)
sample_group <- as.data.frame(c(rep("malignant",2000),rep("normal",2000)))
colnames(sample_group) <- "cluster"
library(ComplexHeatmap)
library(circlize)
col_fun = colorRamp2(c(-2.5, 0, 2.5), c("#00FF00", "#3B3B3B", "#EE0000"))
top_anno <- HeatmapAnnotation(Cluster = sample_group$cluster,
                              col = list(Cluster = c("malignant"= "#F8766D","normal"= "#00BFC4"),border = TRUE))
column_split = sample_group$cluster


ComplexHeatmap::Heatmap(data4,cluster_rows = F,cluster_columns = F,name = " ",
                        show_column_names = F,show_row_names = T,show_heatmap_legend = F,
                        col = col_fun,column_split = column_split,row_title = "Pathway")







a <- fread("./model_built/pathway_score_213/model_result/cell_statistics.csv")
b <- a[5:8,]

data1 <- melt(b)
data1$sample <- factor(data1$sample,levels = c("ATC","TNBC","IDC","DCIS"))
ggplot(data=data1) +
  geom_bar(aes(x=sample, y=value, fill=variable), 
           stat="identity")+
  scale_fill_manual(values = c("#00BFC4","#F8766D"))+
  theme_prism()+
  labs(y="Number of cells")+
  theme(axis.title.x = element_blank())




