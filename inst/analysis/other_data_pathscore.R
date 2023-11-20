
##GSE151530
remove(list = ls())
setwd("~/project/mcIdentify/data/")
library(data.table)
library(dplyr)
infor_data1 <- fread("./model_built/datasets/GSE151530_anno.txt")
expr_data1 <- fread("./model_built/datasets/GSE151530_tpm2.txt")

expr_data <- expr_data1
rownames(expr_data) <- expr_data$V1

expr_matrix <- expr_data[,-1]
expr_matrix <- as.data.frame(t(expr_matrix))
colnames(expr_matrix) <- rownames(expr_data)

data1 <- expr_matrix

GSE256_diff_path3 <- readRDS("./model_built/GSE673_diff_path.rds")
pathway_gene <- readRDS("./KEGG_pathway_gene.rds")

score_gene <- pathway_gene %>% filter(hsa %in% GSE256_diff_path3$hsa)

##pathway score

myFun1 <- function(a){
  
  sum(a)/length(a)
  
}

all_pathway_score <- NA
for (i in 1:213) {
  
  gene <- score_gene %>% filter(hsa == names(table(score_gene$hsa))[i])
  
  a <- data1 %>% select(gene$gene_id[which(gene$gene_id %in% colnames(data1)==TRUE)])
  
  path_score <- as.data.frame(apply(a, 1, myFun1))
  colnames(path_score) <- names(table(score_gene$hsa))[i]
  
  all_pathway_score <- cbind(all_pathway_score,path_score)
  
}

pathway_score <- all_pathway_score
pathway_score <- pathway_score[,-1]

##test
diff_path <- pathway_score
diff_path$Cell <- rownames(diff_path)
infor_data2 <- infor_data1 %>% filter(Type != "unclassified")
diff_path <- left_join(infor_data2,diff_path)
diff_path <- diff_path %>% select(type,GSE256_diff_path3$hsa)

fwrite(diff_path,"./model_built/pathway_score/GSE530_pathway_score.csv")









##GSE146771
remove(list = ls())
setwd("~/project/mcIdentify/data/")
library(data.table)
library(dplyr)
infor_data1 <- fread("./model_built/datasets/GSE146771_anno.txt")
expr_data1 <- fread("./model_built/datasets/GSE146771_tpm.txt")

expr_data <- expr_data1
rownames(expr_data) <- expr_data$V1

expr_matrix <- expr_data[,-1]
expr_matrix <- as.data.frame(t(expr_matrix))
colnames(expr_matrix) <- rownames(expr_data)

data1 <- expr_matrix

GSE256_diff_path3 <- readRDS("./model_built/GSE673_diff_path.rds")
pathway_gene <- readRDS("./KEGG_pathway_gene.rds")

score_gene <- pathway_gene %>% filter(hsa %in% GSE256_diff_path3$hsa)
##pathway score
myFun1 <- function(a){
  
  sum(a)/length(a)
  
}

all_pathway_score <- NA
for (i in 1:213) {
  
  gene <- score_gene %>% filter(hsa == names(table(score_gene$hsa))[i])
  
  a <- data1 %>% select(gene$gene_id[which(gene$gene_id %in% colnames(data1)==TRUE)])
  
  path_score <- as.data.frame(apply(a, 1, myFun1))
  colnames(path_score) <- names(table(score_gene$hsa))[i]
  
  all_pathway_score <- cbind(all_pathway_score,path_score)
  
}

pathway_score <- all_pathway_score
pathway_score <- pathway_score[,-1]


##test
diff_path <- pathway_score
diff_path$CellName <- rownames(diff_path)
diff_path <- left_join(infor_data1,diff_path)
diff_path <- diff_path %>% filter(type == "malignant" | type == "normal") %>%
  select(type,GSE256_diff_path3$hsa)

fwrite(diff_path,"./model_built/pathway_score/GSE771_pathway_score.csv")








remove(list = ls())
##GOSH
setwd("~/project/mcIdentify/data/")
library(data.table)
library(dplyr)
infor_data1 <- fread("./model_built/datasets/GOSH_anno.txt")
expr_data1 <- fread("./model_built/datasets/GOSH_tpm.txt")

##PMC
infor_data1 <- fread("./model_built/datasets/PMC_anno.txt")
expr_data1 <- fread("./model_built/datasets/PMC_tpm.txt")


GSE256_diff_path3 <- readRDS("./model_built/GSE673_diff_path.rds")
pathway_gene <- readRDS("./KEGG_pathway_gene.rds")


expr_data <- expr_data1
rownames(expr_data) <- expr_data1$V1
expr_matrix <- expr_data[,-1]
expr_matrix <- as.data.frame(t(expr_matrix))
colnames(expr_matrix) <- rownames(expr_data) 

data1 <- expr_matrix
score_gene <- pathway_gene %>% filter(hsa %in% GSE256_diff_path3$hsa)


##pathway score

myFun1 <- function(a){
  
  sum(a)/length(a)
  
}

all_pathway_score <- NA
for (i in 1:213) {
  
  gene <- score_gene %>% filter(hsa == names(table(score_gene$hsa))[i])
  
  a <- data1 %>% select(gene$gene_id[which(gene$gene_id %in% colnames(data1)==TRUE)])
  
  path_score <- as.data.frame(apply(a, 1, myFun1))
  colnames(path_score) <- names(table(score_gene$hsa))[i]
  
  all_pathway_score <- cbind(all_pathway_score,path_score)
  
}

pathway_score <- all_pathway_score
pathway_score <- pathway_score[,-1]

##test
diff_path <- pathway_score
diff_path$V1 <- rownames(diff_path)
infor_data2 <- infor_data1 %>% filter(type != "unclassified")
diff_path <- left_join(infor_data2,diff_path)
diff_path <- diff_path %>% select(type,GSE256_diff_path3$hsa)
fwrite(diff_path,"./model_built/pathway_score/GOSH_pathway_score.csv")

fwrite(diff_path,"./model_built/pathway_score/PMC_pathway_score.csv")






## GSE131309 Seq
remove(list = ls())
setwd("~/project/mcIdentify/data/")
library(data.table)
library(dplyr)
infor_data1 <- fread("./model_built/datasets/GSE131309_Seq_anno.txt")
expr_data1 <- fread("./model_built/datasets/GSE131309_Seq_tpm.txt")


GSE256_diff_path3 <- readRDS("./model_built/GSE673_diff_path.rds")
pathway_gene <- readRDS("./KEGG_pathway_gene.rds")

expr_data <- expr_data1
rownames(expr_data) <- expr_data1$V1

expr_matrix <- expr_data[,-1]
expr_matrix <- as.data.frame(t(expr_matrix))
colnames(expr_matrix) <- rownames(expr_data) 
data1 <- expr_matrix
score_gene <- pathway_gene %>% filter(hsa %in% GSE256_diff_path3$hsa)

##pathway score
myFun1 <- function(a){
  
  sum(a)/length(a)
  
}

all_pathway_score <- NA
for (i in 1:213) {
  
  gene <- score_gene %>% filter(hsa == names(table(score_gene$hsa))[i])
  
  a <- data1 %>% select(gene$gene_id[which(gene$gene_id %in% colnames(data1)==TRUE)])
  
  path_score <- as.data.frame(apply(a, 1, myFun1))
  colnames(path_score) <- names(table(score_gene$hsa))[i]
  
  all_pathway_score <- cbind(all_pathway_score,path_score)
  
}

pathway_score <- all_pathway_score
pathway_score <- pathway_score[,-1]

##test
diff_path <- pathway_score
diff_path$barcode <- rownames(diff_path)
infor_data2 <- infor_data1 %>% filter(type != "unclassified")
diff_path <- left_join(infor_data2,diff_path)
diff_path <- diff_path %>% select(type,GSE256_diff_path3$hsa)

fwrite(diff_path,"./model_built/pathway_score/GSE309_Seq_pathway_score.csv")
fwrite(diff_path,"./model_built/pathway_score/GSE309_10X_pathway_score.csv")





