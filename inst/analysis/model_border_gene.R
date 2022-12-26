

setwd("~/project/mcIdentify/data/")
remove(list = ls())
library(data.table)
library(dplyr)

data1 <- fread("./processed_data/GSE151530_tpm.txt")
expr_data <- as.data.frame(data1[,-1])
rownames(expr_data) <- data1$V1



## pathway score
setwd("~/project/mcIdentify/data/")
remove(list = ls())
library(data.table)
library(dplyr)

data <- fread("./model_built/pathway_score_213/border_data_gene/GSE151530_tpm_500.txt")
data <- fread("./model_built/pathway_score_213/border_data_gene/GSE151530_tpm_700.txt")


infor_data1 <- fread("./model_built/datasets/GSE151530_anno.txt")
pathway_gene <- readRDS("./KEGG_pathway_gene.rds")
expr_data1 <- data %>% filter(V1 %in% names(table(pathway_gene$gene_id)))


expr_data <- expr_data1
rownames(expr_data) <- expr_data$V1

expr_matrix <- expr_data[,-1]
expr_matrix <- as.data.frame(t(expr_matrix))
colnames(expr_matrix) <- rownames(expr_data)

data1 <- expr_matrix

GSE256_diff_path3 <- readRDS("./model_built/pathway_score_213/GSE673_diff_path213.rds")

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


diff_path <- pathway_score
diff_path$Cell <- rownames(diff_path)
infor_data2 <- infor_data1 %>% filter(Type != "unclassified")
diff_path <- left_join(diff_path,infor_data2)
diff_path <- diff_path %>% select(type,GSE256_diff_path3$hsa) %>% na.omit()


fwrite(diff_path,"./model_built/pathway_score_213/border_data_gene/pathway_score/GSE530_500.csv")
fwrite(diff_path,"./model_built/pathway_score_213/border_data_gene/pathway_score/GSE530_700.csv")





for (number in c(900,1100,1300,1500)) {
  read_filename <- paste0("./model_built/pathway_score_213/border_data_gene/GSE151530_tpm_",number,".txt")
  data <- fread(read_filename)

  
  infor_data1 <- fread("./model_built/datasets/GSE151530_anno.txt")
  pathway_gene <- readRDS("./KEGG_pathway_gene.rds")
  expr_data1 <- data %>% filter(V1 %in% names(table(pathway_gene$gene_id)))
  
  
  expr_data <- expr_data1
  rownames(expr_data) <- expr_data$V1
  
  expr_matrix <- expr_data[,-1]
  expr_matrix <- as.data.frame(t(expr_matrix))
  colnames(expr_matrix) <- rownames(expr_data)
  data1 <- expr_matrix
  
  GSE256_diff_path3 <- readRDS("./model_built/pathway_score_213/GSE673_diff_path213.rds")
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
  

  diff_path <- pathway_score
  diff_path$Cell <- rownames(diff_path)
  infor_data2 <- infor_data1 %>% filter(Type != "unclassified")
  diff_path <- left_join(diff_path,infor_data2)
  diff_path <- diff_path %>% select(type,GSE256_diff_path3$hsa) %>% na.omit()
  
  
  write_filename <- paste0("./model_built/pathway_score_213/border_data_gene/pathway_score/GSE530_",number,".csv")
  fwrite(diff_path,write_filename)
}




###gene+pathway
setwd("~/project/mcIdentify/data/")
remove(list = ls())
library(data.table)
library(dplyr)

all_data <- fread("./processed_data/GSE148673_tpm.txt")
border_data <- as.data.frame(all_data[,-1])
rownames(border_data) <- all_data$V1


for (number in c(500,1000,1500,2000,2500)) {

##border gene select
low_number <- NA
testdata <- border_data
for (i in 1:ncol(border_data)) {
  gene_number <- length(which(testdata[,i] > 0))
  judge <- gene_number - number
  
  if (judge >= 0) {
    random_number <- sample(1:gene_number, judge, replace = FALSE)
    testdata[which(testdata[,i] > 0)[random_number],i] <- 0
  }else{
    low_number <- append(low_number,i)
  }
  
}

if (number == 500) {
  testdata1 <- testdata
}else{
  testdata1 <- testdata[,-low_number[-1]] 
}


write_filename_border <- paste0("./model_built/pathway_score_213/border_data_gene/GSE148673_tpm_",number,".txt")

fwrite(testdata1,write_filename_border,row.names = T)



####pathway score
read_filename <- paste0("./model_built/pathway_score_213/border_data_gene/GSE148673_tpm_",number,".txt")
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

GSE256_diff_path3 <- readRDS("./model_built/pathway_score_213/GSE673_diff_path213.rds")
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


diff_path <- pathway_score
diff_path$Cell <- rownames(diff_path)
# infor_data2 <- infor_data1 %>% filter(Type != "unclassified")
infor_data2 <- infor_data1
diff_path <- left_join(diff_path,infor_data2)
diff_path <- diff_path %>% select(type,GSE256_diff_path3$hsa) %>% na.omit()


write_filename <- paste0("./model_built/pathway_score_213/border_data_gene/pathway_score/GSE673/GSE673_",number,".csv")
fwrite(diff_path,write_filename)
}





##GSE530
a <- NA
for (i in 1:ncol(border_data)) {
  
  a <- append(a,table(border_data[,i]>0))
  
}

a <- as.data.frame(a[-1])
b <- as.data.frame(a[seq(2,nrow(a),2),])

colnames(b) <- "gene_number"

border_500_1000 <- border_data[,which(b$gene_number>=500 & b$gene_number <1000)]
border_1000_1500 <- border_data[,which(b$gene_number>=1000 & b$gene_number <1500)]
border_1500_2000 <- border_data[,which(b$gene_number>=1500 & b$gene_number <2000)]
border_2000_2500 <- border_data[,which(b$gene_number>=2000 & b$gene_number <2500)]
border_2500_00 <- border_data[,which(b$gene_number>=2500)]


fwrite(border_500_1000,"./model_built/pathway_score_213/border_interval/border_500_1000.txt",row.names = T)
fwrite(border_1000_1500,"./model_built/pathway_score_213/border_interval/border_1000_1500.txt",row.names = T)
fwrite(border_1500_2000,"./model_built/pathway_score_213/border_interval/border_1500_2000.txt",row.names = T)
fwrite(border_2000_2500,"./model_built/pathway_score_213/border_interval/border_2000_2500.txt",row.names = T)
fwrite(border_2500_00,"./model_built/pathway_score_213/border_interval/border_2500_00.txt",row.names = T)



file_name <- list.files("./model_built/pathway_score_213/border_interval/")[1:5]
for (number in file_name) {
  
  read_filename <- paste0("./model_built/pathway_score_213/border_interval/",number)
  data <- fread(read_filename)
  
  
  infor_data1 <- fread("./model_built/datasets/GSE151530_anno.txt")
  pathway_gene <- readRDS("./KEGG_pathway_gene.rds")
  expr_data1 <- data %>% filter(V1 %in% names(table(pathway_gene$gene_id)))
  
  
  expr_data <- expr_data1
  rownames(expr_data) <- expr_data$V1
  
  expr_matrix <- expr_data[,-1]
  expr_matrix <- as.data.frame(t(expr_matrix))
  colnames(expr_matrix) <- rownames(expr_data)
  data1 <- expr_matrix
  
  GSE256_diff_path3 <- readRDS("./model_built/pathway_score_213/GSE673_diff_path213.rds")
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
  

  diff_path <- pathway_score
  diff_path$Cell <- rownames(diff_path)
  infor_data2 <- infor_data1 %>% filter(Type != "unclassified")
  diff_path <- left_join(diff_path,infor_data2)
  diff_path <- diff_path %>% select(type,GSE256_diff_path3$hsa) %>% na.omit()
  
  
  write_filename <- paste0("./model_built/pathway_score_213/border_interval/pathway_score_interval/",number,".csv")
  fwrite(diff_path,write_filename)
}





setwd("~/project/mcIdentify/data/")
remove(list = ls())
library(data.table)
library(dplyr)


all_data <- fread("./model_built/pathway_score_213/border_interval/border_2500_00.txt")
border_data <- as.data.frame(all_data[,-1])
rownames(border_data) <- all_data$V1



for (number in c(500,1000,1500,2000,2500)) {
  
  ##border gene select
  low_number <- NA
  testdata <- border_data
  for (i in 1:ncol(border_data)) {
    gene_number <- length(which(testdata[,i] > 0))
    judge <- gene_number - number
    
    if (judge >= 0) {
      random_number <- sample(1:gene_number, judge, replace = FALSE)
      testdata[which(testdata[,i] > 0)[random_number],i] <- 0
    }
    
  }
  testdata1 <- testdata
  
  write_filename_border <- paste0("./model_built/pathway_score_213/border_data_gene/GSE151530_tpm_",number,".txt")
  
  fwrite(testdata1,write_filename_border,row.names = T)
  
  
  
  ####pathway score
  read_filename <- paste0("./model_built/pathway_score_213/border_data_gene/GSE151530_tpm_",number,".txt")
  data <- fread(read_filename)
  
  
  infor_data1 <- fread("./model_built/datasets/GSE151530_anno.txt")
  pathway_gene <- readRDS("./KEGG_pathway_gene.rds")
  expr_data1 <- data %>% filter(V1 %in% names(table(pathway_gene$gene_id)))
  
  
  expr_data <- expr_data1
  rownames(expr_data) <- expr_data$V1
  
  expr_matrix <- expr_data[,-1]
  expr_matrix <- as.data.frame(t(expr_matrix))
  colnames(expr_matrix) <- rownames(expr_data)
  data1 <- expr_matrix
  
  GSE256_diff_path3 <- readRDS("./model_built/pathway_score_213/GSE673_diff_path213.rds")
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
  

  diff_path <- pathway_score
  diff_path$Cell <- rownames(diff_path)
  infor_data2 <- infor_data1 %>% filter(Type != "unclassified")
  diff_path <- left_join(diff_path,infor_data2)
  diff_path <- diff_path %>% select(type,GSE256_diff_path3$hsa) %>% na.omit()
  
  
  write_filename <- paste0("./model_built/pathway_score_213/border_data_gene/pathway_score/GSE530_",number,".csv")
  fwrite(diff_path,write_filename)
}
 



#### other sample GOSH

setwd("~/project/mcIdentify/data/")
remove(list = ls())
library(data.table)
library(dplyr)

all_data <- fread("./model_built/datasets/GOSH_tpm.txt")
border_data <- as.data.frame(all_data[,-1])
rownames(border_data) <- all_data$V1


a <- NA
for (i in 1:ncol(border_data)) {
  
  a <- append(a,table(border_data[,i]>0))
  
}

a <- as.data.frame(a[-1])
b <- as.data.frame(a[seq(2,nrow(a),2),])

colnames(b) <- "gene_number"

border_500 <- border_data[,which(b$gene_number < 500)]
border_500_1000 <- border_data[,which(b$gene_number>=500 & b$gene_number <1000)]
border_1000_1500 <- border_data[,which(b$gene_number>=1000 & b$gene_number <1500)]
border_1500_2000 <- border_data[,which(b$gene_number>=1500 & b$gene_number <2000)]
border_2000_2500 <- border_data[,which(b$gene_number>=2000 & b$gene_number <2500)]
border_2500_00 <- border_data[,which(b$gene_number>=2500)]


fwrite(border_500,"./model_built/pathway_score_213/border_interval/GOSH/G0SH_border_500.txt",row.names = T)
fwrite(border_500_1000,"./model_built/pathway_score_213/border_interval/GOSH/G0SH_border_500_1000.txt",row.names = T)
fwrite(border_1000_1500,"./model_built/pathway_score_213/border_interval/GOSH/G0SH_border_1000_1500.txt",row.names = T)
fwrite(border_1500_2000,"./model_built/pathway_score_213/border_interval/GOSH/G0SH_border_1500_2000.txt",row.names = T)
fwrite(border_2000_2500,"./model_built/pathway_score_213/border_interval/GOSH/G0SH_border_2000_2500.txt",row.names = T)
fwrite(border_2500_00,"./model_built/pathway_score_213/border_interval/GOSH/G0SH_border_2500_00.txt",row.names = T)


file_name <- list.files("./model_built/pathway_score_213/border_interval/GOSH/")
for (number in file_name) {
  
  read_filename <- paste0("./model_built/pathway_score_213/border_interval/GOSH/",number)
  data <- fread(read_filename)
  
  
  infor_data1 <- fread("./model_built/datasets/GOSH_anno.txt")
  pathway_gene <- readRDS("./KEGG_pathway_gene.rds")
  expr_data1 <- data %>% filter(V1 %in% names(table(pathway_gene$gene_id)))
  
  
  expr_data <- expr_data1
  rownames(expr_data) <- expr_data$V1
  
  expr_matrix <- expr_data[,-1]
  expr_matrix <- as.data.frame(t(expr_matrix))
  colnames(expr_matrix) <- rownames(expr_data)
  data1 <- expr_matrix
  
  GSE256_diff_path3 <- readRDS("./model_built/pathway_score_213/GSE673_diff_path213.rds")
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
  

  diff_path <- pathway_score
  diff_path$V1 <- rownames(diff_path)
  infor_data2 <- infor_data1 %>% filter(type != "unclassified")
  diff_path <- left_join(diff_path,infor_data2)
  diff_path <- diff_path %>% select(type,GSE256_diff_path3$hsa) %>% na.omit()
  
  
  write_filename <- paste0("./model_built/pathway_score_213/border_interval/pathway_score_interval/",number,".csv")
  fwrite(diff_path,write_filename)
}








#### other sample GOSH

setwd("~/project/mcIdentify/data/")
remove(list = ls())
library(data.table)
library(dplyr)

all_data <- fread("./model_built/datasets/GSE148673_tpm.txt")
border_data <- as.data.frame(all_data[,-1])
rownames(border_data) <- all_data$V1


a <- NA
for (i in 1:ncol(border_data)) {
  
  a <- append(a,table(border_data[,i]>0))
  
}

a <- as.data.frame(a[-1])
b <- as.data.frame(a[seq(2,nrow(a),2),])

colnames(b) <- "gene_number"
border_500 <- border_data[,which(b$gene_number < 500)]
border_500_1000 <- border_data[,which(b$gene_number>=500 & b$gene_number <1000)]
border_1000_1500 <- border_data[,which(b$gene_number>=1000 & b$gene_number <1500)]
border_1500_2000 <- border_data[,which(b$gene_number>=1500 & b$gene_number <2000)]
border_2000_2500 <- border_data[,which(b$gene_number>=2000 & b$gene_number <2500)]
border_2500_00 <- border_data[,which(b$gene_number>=2500)]


fwrite(border_500,"./model_built/pathway_score_213/border_interval/GSE673/GSE673_border_500.txt",row.names = T)
fwrite(border_500_1000,"./model_built/pathway_score_213/border_interval/GSE673/GSE673_border_500_1000.txt",row.names = T)
fwrite(border_1000_1500,"./model_built/pathway_score_213/border_interval/GSE673/GSE673_border_1000_1500.txt",row.names = T)
fwrite(border_1500_2000,"./model_built/pathway_score_213/border_interval/GSE673/GSE673_border_1500_2000.txt",row.names = T)
fwrite(border_2000_2500,"./model_built/pathway_score_213/border_interval/GSE673/GSE673_border_2000_2500.txt",row.names = T)
fwrite(border_2500_00,"./model_built/pathway_score_213/border_interval/GSE673/GSE673_border_2500_00.txt",row.names = T)


file_name <- list.files("./model_built/pathway_score_213/border_interval/GSE673/")
for (number in file_name) {
  
  read_filename <- paste0("./model_built/pathway_score_213/border_interval/GSE673/",number)
  data <- fread(read_filename)
  
  
  infor_data1 <- fread("./model_built/datasets/GSE148673_anno.txt")
  pathway_gene <- readRDS("./KEGG_pathway_gene.rds")
  expr_data1 <- data %>% filter(V1 %in% names(table(pathway_gene$gene_id)))
  
  
  expr_data <- expr_data1
  rownames(expr_data) <- expr_data$V1
  
  expr_matrix <- expr_data[,-1]
  expr_matrix <- as.data.frame(t(expr_matrix))
  colnames(expr_matrix) <- rownames(expr_data)
  data1 <- expr_matrix
  
  GSE256_diff_path3 <- readRDS("./model_built/pathway_score_213/GSE673_diff_path213.rds")
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

  diff_path <- pathway_score
  diff_path$barcode <- rownames(diff_path)
  infor_data2 <- infor_data1 %>% filter(type != "unclassified")
  diff_path <- left_join(diff_path,infor_data2)
  diff_path <- diff_path %>% select(type,GSE256_diff_path3$hsa) %>% na.omit()
  
  
  write_filename <- paste0("./model_built/pathway_score_213/border_interval/pathway_score_interval/",number,".csv")
  fwrite(diff_path,write_filename)
}




