

setwd("~/project/mcIdentify/data/")


library(data.table)
library(dplyr)

path_impor1 <- fread("./model_built/pathway_score_213/pathway_importance/path_impor_GSE673_model15.csv")
path_score1 <- fread("./model_built/pathway_score_213/GSE673_pathway_score213.csv")
path_impor1$hsa <- colnames(path_score1)

path_impor1 <- path_impor1[-1,-1]
path_impor1$number <- c(1:213)

plot(path_impor1$V2,ylim = c(0.048,0.11),col = "blue", pch = 19, cex = 1)


path_impor1 <- path_impor1 %>% dplyr::mutate(fac = case_when(V2 > 0.058 ~ "A", TRUE ~ "B"))
path_impor1$fac <- as.factor(path_impor1$fac)

path_impor1 <- left_join(path_impor1,GSE673_diff_path213)

gene_list <- KEGG_pathway_gene %>% filter(hsa %in% path_impor1[path_impor1$fac == "A",]$hsa)


fwrite(gene_list,"./model_built/pathway_score_213/pathway_importance/gene_list.csv")



x <- list("Oxidative phosphorylation" = gene_list[gene_list$hsa=="hsa00190",]$gene_id,
          
          "Viral myocarditis" = gene_list[gene_list$hsa=="hsa05416",]$gene_id,
          
          "Type I diabetes mellitus" = gene_list[gene_list$hsa=="hsa04940",]$gene_id,
          
          "Antigen processing and presentation" = gene_list[gene_list$hsa=="hsa04612",]$gene_id)


venn.plot <- venn.diagram(
  x,
  filename = NULL,
  lty = 1,
  lwd = 1,
  col = "black",  
  fill = c("#6E9ECE", "#EFDBB9","#E6928F","4E9595"),
  alpha = 0.60,
  cat.col = "black",
  cat.cex = 0.8,
  cat.fontface = "bold",
  margin = 0.07,
  cex = 0.8
)




pdf("venn.pdf",width = 12,height = 12,pointsize = 20.5)
grid.draw(venn.plot)
dev.off()





gene1 <- as.data.frame(sort(table(gene_list$gene_id)))
gene2 <- gene1[gene1$Freq > 2,]


data1 <- fread("./model_built/datasets/GSE148673_tpm.txt")
data2 <- data1 %>% filter(V1 %in% gene2$Var1)

data3 <- data2[,-1]
data4 <- as.data.frame(t(data3))
colnames(data4) <- data2$V1
data4$barcode <- rownames(data4)

infor <- fread("./model_built/datasets/GSE148673_anno.txt")
infor <- infor %>% mutate(type= if_else(cluster.pred == "T","malignant","normal"))

infor1 <- infor %>% select(barcode,type)

data5 <- left_join(infor1,data4)
data5 <- data5 %>% select(-barcode)

data6 <- melt(data5)

data6$variable <- factor(data6$variable,levels = c("HLA-A","HLA-B","HLA-C","HLA-E","HLA-F","HLA-G",
                                                   "HLA-DRA","HLA-DRB1","HLA-DRB5","HLA-DQA1","HLA-DQA2",
                                                   "HLA-DQB1","HLA-DOB","HLA-DMA","HLA-DMB","HLA-DOA","HLA-DPA1","HLA-DPB1"))

library(ggpubr)
library(ggprism)
library(ggplot2)
library(cowplot)



ggplot(data=data6,aes(x=variable,y=value,fill=factor(type)))+
  geom_violin()+
  stat_compare_means(aes(label = ..p.signif..))+
  theme_prism()+
  theme(axis.text.x = element_text(angle = 30,vjust = 1, hjust = 1) )+
  labs(y="Gene expression",title = "")+
  theme(axis.title.x = element_blank())


