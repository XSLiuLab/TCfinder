

setwd("~/project/mcIdentify/data/")
remove(list = ls())
library(data.table)
library(dplyr)

all_data <- fread("./processed_data/GSE148673_tpm.txt")
border_data <- as.data.frame(all_data[,-1])
rownames(border_data) <- all_data$V1




a <- apply(border_data, 2, function(x){which(x > 0)})

for (number in c(20)) {
  
  ##border gene select
  low_number <- NA
  testdata <- border_data
  for (i in 1:35727) {
    random_number <- sample(a[[i]], number, replace = FALSE)
    testdata[random_number,i] <- 0

  }
  



  write_filename_border <- paste0("./model_built/pathway_score_213/simulation_gene/GSE148673_tpm_",number,".txt")
  
  fwrite(testdata,write_filename_border,row.names = T)
  
  
  
  ####pathway score
  read_filename <- paste0("./model_built/pathway_score_213/simulation_gene/GSE148673_tpm_",number,".txt")
  data <- fread(read_filename)
  
  
  infor_data1 <- fread("./model_built/datasets/GSE148673_anno.txt")
  infor_data1 <- infor_data1 %>% mutate(type = case_when(cluster.pred == "T"~"malignant",
                                                         cluster.pred == "N"~"normal"))
  
  pathway_gene <- readRDS("./KEGG_pathway_gene.rds")
  expr_data1 <- data %>% filter(V1 %in% names(table(pathway_gene$gene_id)))
  
  
  expr_data <- expr_data1
  rownames(expr_data) <- expr_data$V1
  
  expr_matrix <- expr_data[,-1]
  expr_matrix <- as.data.frame(t(expr_matrix))
  colnames(expr_matrix) <- rownames(expr_data)
  data1 <- expr_matrix
  
  GSE673_diff_path3 <- readRDS("./model_built/pathway_score_213/GSE673_diff_path213.rds")
  score_gene <- pathway_gene %>% filter(hsa %in% GSE673_diff_path3$hsa)
  
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

  diff_path <- pathway_score
  diff_path$barcode <- rownames(diff_path)
  # infor_data2 <- infor_data1 %>% filter(Type != "unclassified")
  infor_data2 <- infor_data1
  diff_path <- left_join(diff_path,infor_data2)
  diff_path <- diff_path %>% select(type,GSE256_diff_path3$hsa) %>% na.omit()
  
  
  write_filename <- paste0("./model_built/pathway_score_213/simulation_gene/pathway_score/GSE673_",number,".csv")
  fwrite(diff_path,write_filename)
}




