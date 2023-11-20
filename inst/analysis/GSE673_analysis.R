

setwd("~/project/mcIdentify/data/")
remove(list = ls())
library(data.table)
library(dplyr)

GSE673_eff <- fread("./model_built/pathway_score_213/model_result/GSE673_inter_eff.csv")

tumor_type <- GSE673_eff[1:4,]

tumor_type1 <- melt(tumor_type)

V1 <- c("GSE673_ATC","GSE673_IDC","GSE673_TNBC","GSE673_DCIS")
tumor_type1$V1 <- factor(tumor_type1$V1,levels = c("GSE673_ATC","GSE673_IDC","GSE673_TNBC","GSE673_DCIS"))
library(ggplot2)
library(ggprism)
ggplot(tumor_type1,aes(x=V1,y=value,fill=variable))+
  geom_bar(position="dodge",stat="identity")+
  labs(x="Cancer type",y="")+
  theme_prism()+
  geom_hline(aes(yintercept=0.95),linetype=5,col="black")+
  scale_y_continuous(breaks=c(0,.25,0.5,0.75,0.95,1))+
  scale_x_discrete(breaks=V1, labels=c("ATC","IDC","TNBC","DCIS"))+
  scale_fill_manual(values = c("#6E9ECE", "#CCCCCC","#E6928F","#8FBC8F"),
                    breaks=c("f1", "accuracy", "recall", "precisoon"),
                    labels=c("F1 score", "Accuracy", "Recall", "Precison"))




###gene number
gene_number <- GSE673_eff[5:10,]
gene_number1 <- melt(gene_number)
V1 <- c("GSE673_500","GSE673_500_1000","GSE673_1000_1500","GSE673_1500_2000","GSE673_2000_2500","GSE673_2500_00")
gene_number1$V1 <- factor(gene_number1$V1,levels = V1)
library(ggplot2)
library(ggprism)
ggplot(gene_number1,aes(x=V1,y=value,fill=variable))+
  geom_bar(position="dodge",stat="identity")+
  labs(x="Gene number",y="")+
  theme_prism()+
  geom_hline(aes(yintercept=0.95),linetype=5,col="black")+
  scale_y_continuous(breaks=c(0,.25,0.5,0.75,0.95,1))+
  scale_x_discrete(breaks=V1, labels=c("<500","500~1000","1000~1500","1500~2000","2000~2500",">2500"))+
  scale_fill_manual(values = c("#6E9ECE", "#CCCCCC","#E6928F","#8FBC8F"),
                    breaks=c("f1", "accuracy", "recall", "precisoon"),
                    labels=c("F1 score", "Accuracy", "Recall", "Precison"))



# simulation gene
gene_number <- GSE673_eff[11:15,]
gene_number1 <- melt(gene_number)
V1 <- c("simulation_500","simulation_1000","simulation_1500","simulation_2000","simulation_2500")
gene_number1$V1 <- factor(gene_number1$V1,levels = V1)
library(ggplot2)
library(ggprism)
ggplot(gene_number1,aes(x=V1,y=value,fill=variable))+
  geom_bar(position="dodge",stat="identity")+
  labs(x="Simulate gene number",y="")+
  theme_prism()+
  geom_hline(aes(yintercept=0.95),linetype=5,col="black")+
  scale_y_continuous(breaks=c(0,.25,0.5,0.75,0.95,1))+
  scale_x_discrete(breaks=V1, labels=c("500","1000","1500","2000","2500"))+
  scale_fill_manual(values = c("#6E9ECE", "#CCCCCC","#E6928F","#8FBC8F"),
                    breaks=c("f1", "accuracy", "recall", "precisoon"),
                    labels=c("F1 score", "Accuracy", "Recall", "Precison"))




# setwd("~/project/mcIdentify/data/")
# remove(list = ls())
# 
# library(data.table)
# library(dplyr)
# GSE673_diff_path213 <- readRDS("~/project/mcIdentify/data/model_built/pathway_score_213/GSE673_diff_path213.rds")
# data1 <- fread("./model_built/pathway_score_213/tumor_type_data/GSE673_tumor_type_data.csv")
# 
# data2 <- data1 %>% mutate(cancer_type = case_when(cell_type == "ATC1." ~ "ATC",
#                                                   cell_type == "ATC2." ~ "ATC",
#                                                   cell_type == "ATC3." ~ "ATC",
#                                                   cell_type == "ATC4." ~ "ATC",
#                                                   cell_type == "ATC5." ~ "ATC",
#                                                   cell_type == "DCIS1" ~ "DCIS",
#                                                   cell_type == "IDC1." ~ "IDC",
#                                                   cell_type == "IDC2." ~ "IDC",
#                                                   cell_type == "TNBC1" ~ "TNBC",
#                                                   cell_type == "TNBC2" ~ "TNBC",
#                                                   cell_type == "TNBC3" ~ "TNBC",
#                                                   cell_type == "TNBC4" ~ "TNBC",
#                                                   cell_type == "TNBC5" ~ "TNBC"))
# 
# data3 <- data2 %>% select(type,cancer_type,GSE673_diff_path213$hsa)
# single_cancer <- data3 %>% filter(cancer_type == "TNBC") %>% select(type,GSE673_diff_path213$hsa)
# fwrite(single_cancer,"./model_built/pathway_score_213/tumor_type_data/GSE673_TNBC.csv")

