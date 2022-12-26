
setwd("~/project/mcIdentify/data/")
remove(list = ls())

library(data.table)
library(dplyr)


tcga_counts <- readRDS("~/project/common_data/tcga/tcga_clean_counts.rds")

data1 <- tcga_counts %>% apply(2,function(x){x/sum(x) * 10000})

data1 <- as.data.frame(data1)
data2 <- data1 %>% mutate_all(funs(log2(.+1)))
data2 <- round(data2,3)

data3 <- as.data.frame(t(data2))
colnames(data3) <- rownames(tcga_counts)


# pathway socre

GSE256_diff_path3 <- readRDS("./model_built/pathway_score_213/GSE673_diff_path213.rds")
pathway_gene <- readRDS("./KEGG_pathway_gene.rds")
score_gene <- pathway_gene %>% filter(hsa %in% GSE256_diff_path3$hsa)


myFun1 <- function(a){
  
  sum(a)/length(a)
  
}

all_pathway_score <- NA
for (i in 1:213) {
  
  gene <- score_gene %>% filter(hsa == names(table(score_gene$hsa))[i])
  
  a <- data3 %>% select(gene$gene_id[which(gene$gene_id %in% colnames(data3)==TRUE)])
  
  path_score <- as.data.frame(apply(a, 1, myFun1))
  colnames(path_score) <- names(table(score_gene$hsa))[i]
  
  all_pathway_score <- cbind(all_pathway_score,path_score)
  
}

pathway_score <- all_pathway_score
pathway_score <- pathway_score[,-1]

fwrite(pathway_score,"./model_built/pathway_score_213/bulk_pathway_score/bulk_pathway_score.txt",row.names = T)






###figure
setwd("~/project/mcIdentify/data/")
remove(list = ls())

library(data.table)
library(dplyr)
tcga_score <- fread("./model_built/pathway_score_213/bulk_pathway_score/bulk_pathway_score.txt")

data1 <- tcga_score %>% select(V1,hsa00190,hsa04612,hsa04940,hsa05416,hsa04110)
colnames(data1)[1] <- "tcga_id"

library(NeoEnrichment)
data1$cancer_type <- get_cancer_type(data1$tcga_id)

data2 <- data1 %>% dplyr::mutate(tissue = case_when(grepl("*01$",tcga_id) ~ "tumor",
                                     grepl("*11$",tcga_id) ~ "normal")) %>% na.omit()



normal_data <- data2 %>% filter(tissue == "normal")
normal_name <- as.data.frame(sort(table(normal_data$cancer_type)))
normal_name <- normal_name %>% filter(Freq > 20) %>% arrange(-Freq) %>% filter(Var1 != c("PRAD","KICH"))
rownames(normal_name) <- normal_name$Var1

data3 <- data2 %>% filter(cancer_type %in% normal_name$Var1)

data4 <- melt(data3)

pathway1 <- data4 %>% filter(variable == "hsa05416")

all_sample <- pathway1
all_sample$cancer_type <- "Pan-cancer"
pathway2 <- rbind(pathway1,all_sample)
pathway2$tissue <- factor(pathway2$tissue,levels = c("tumor","normal"))
pathway2$cancer_type <- factor(pathway2$cancer_type,levels = c("Pan-cancer",rownames(normal_name)))



library(ggpubr)
library(ggprism)
library(ggplot2)
library(cowplot)
ggplot(data=pathway2,aes(x=cancer_type,y=value,fill=factor(tissue)))+
  geom_boxplot()+
  stat_compare_means(aes(label = ..p.signif..))+
  #ylim(0.3,1.6)+
  theme_prism()+
  labs(y="Pathway score",title = "hsa05416")+
  theme(axis.title.x = element_blank())+
  theme(axis.text.x = element_text(angle = 15))






path_score2 <- fread("./model_built/pathway_score_213/GSE673_pathway_score213.csv")
infor1 <- fread("./model_built/datasets/GSE148673_anno.txt")

data1 <- cbind(infor1,path_score2)

data2 <- data1 %>% select(type,cell_type,hsa00190,hsa04612,hsa04940,hsa05416,hsa04110)

data2 <- data2 %>% mutate(cancer_type = case_when(cell_type == "ATC1." ~ "ATC",
                                                  cell_type == "ATC2." ~ "ATC",
                                                  cell_type == "ATC3." ~ "ATC",
                                                  cell_type == "ATC4." ~ "ATC",
                                                  cell_type == "ATC5." ~ "ATC",
                                                  cell_type == "DCIS1" ~ "DCIS",
                                                  cell_type == "IDC1." ~ "IDC",
                                                  cell_type == "IDC2." ~ "IDC",
                                                  cell_type == "TNBC1" ~ "TNBC",
                                                  cell_type == "TNBC2" ~ "TNBC",
                                                  cell_type == "TNBC3" ~ "TNBC",
                                                  cell_type == "TNBC4" ~ "TNBC",
                                                  cell_type == "TNBC5" ~ "TNBC"))
data2 <- data2 %>% select(-cell_type)


data3 <- melt(data2)

data4 <- data3 %>% filter(variable == "hsa04110")

library(ggpubr)
library(ggprism)
library(ggplot2)
library(cowplot)
ggplot(data=data4,aes(x=cancer_type,y=value,fill=factor(type)))+
  geom_boxplot()+
  stat_compare_means(aes(label = ..p.signif..))+
  #ylim(0.3,1.6)+
  theme_prism()+
  labs(y="Pathway score",title = "hsa04110")+
  theme(axis.title.x = element_blank())+
  theme(axis.text.x = element_text(angle = 15))





