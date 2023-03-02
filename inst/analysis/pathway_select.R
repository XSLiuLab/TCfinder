



remove(list = ls())
setwd("~/project/mcIdentify/data/")
# model
##GSE148673

library(data.table)
library(dplyr)
infor_data1 <- fread("./model_built/datasets/GSE148673_anno.txt")
expr_data1 <- fread("./model_built/datasets/GSE148673_tpm.txt")


expr_data <- expr_data1
rownames(expr_data) <- expr_data1$V1

expr_matrix <- expr_data[,-1]
expr_matrix <- as.data.frame(t(expr_matrix))
colnames(expr_matrix) <- rownames(expr_data)


data1 <- expr_matrix

##pathway score
pathway_gene <- readRDS("./KEGG_pathway_gene.rds")

infor_data1 <- infor_data1 %>% mutate_at(.vars = "cluster.pred",.funs = funs(ifelse(.=="T","malignant","normal")))



myFun1 <- function(a){
  
  sum(a)/length(a)
  
}


all_pathway_score <- NA
for (i in 1:335) {
  
  gene <- pathway_gene %>% filter(hsa == names(table(pathway_gene$hsa))[i])
  
  a <- data1 %>% select(gene$gene_id[which(gene$gene_id %in% colnames(data1)==TRUE)])
  
  path_score <- as.data.frame(apply(a, 1, myFun1))
  colnames(path_score) <- names(table(pathway_gene$hsa))[i]
  
  all_pathway_score <- cbind(all_pathway_score,path_score)
  
}


pathway_score <- all_pathway_score
pathway_score <- pathway_score[,-1]



##test

diff_path <- pathway_score
diff_path$barcode <- rownames(diff_path)
diff_path <- left_join(infor_data1,diff_path)
diff_path <- diff_path %>% filter(cluster.pred == "malignant" | cluster.pred == "normal") %>%
  select(cluster.pred,names(table(pathway_gene$hsa)))

colnames(diff_path)[1] <- "type"

saveRDS(diff_path,"./model_built/GSE673_all_pathway_score.rds")


###boxplot
library(ggpubr)
library(ggprism)
library(ggplot2)
plot1 <- ggplot(data=diff_path,aes(x=PredictionRefined,y=hsa00010))+
  geom_boxplot(size=1)+
  stat_compare_means(label.x=1.2,size=5,method = "wilcox.test")+
  theme_prism()+
  labs(y="Pathway Score")+
  theme(axis.title.x = element_blank())

plot1


##testing
tumor_sample <- diff_path %>% filter(type=="malignant")
normal_sample <- diff_path %>% filter(type=="normal")

tumor_infor <- infor_data1 %>% filter(cluster.pred=="malignant")


diff_path2 <- diff_path %>% select(-type)

test_pathway <- NA

for (name in colnames(diff_path2)) {
  
  if (sum(select(diff_path2,name) == 0) < nrow(diff_path2)*0.01) {
    
    a <- wilcox.test(as.matrix(select(tumor_sample,name)),as.matrix(select(normal_sample,name)),paired = F, correct = F)
    b <- as.data.frame(a$p.value)
    rownames(b) <- name
    test_pathway <- rbind(test_pathway,b)
    
  }
  
}

test_value <- test_pathway
##select pathway
colnames(test_value) <- "pvalue"
test_value2 <- test_value %>% filter(pvalue< 0.05)
test_value2$hsa <- rownames(test_value2)
test_value2 <- test_value2 %>% arrange(pvalue)



##pathway
path_test <- test_value2
pathway <- pathway_gene %>% filter(!duplicated(hsa)) %>% select(hsa,pathway_id)

path_test <- left_join(path_test,pathway) %>% filter(pvalue == 0)
saveRDS(path_test,"./model_built/GSE673_diff_path.rds")


testdata <- diff_path %>% select(type,path_test$hsa)
fwrite(testdata,"./model_built/GSE673_pathway_score.csv")


