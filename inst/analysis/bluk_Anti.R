

remove(list = ls())
setwd("~/project/mcIdentify/data/")

tcga_counts <- readRDS("~/project/common_data/tcga/tcga_clean_counts.rds")

library(NeoEnrichment)
library(dplyr)
library(data.table)
infor <- as.data.frame(colnames(tcga_counts))
colnames(infor) <- "barcode"
infor$type <- get_cancer_type(infor$barcode)
  
BRCA_infor <- infor %>% filter(type == "BRCA") 
BRCA_infor <- BRCA_infor %>% dplyr::mutate(tissue = case_when(grepl("*01$",BRCA_infor$barcode) ~ "tumor",
                                  grepl("*11$",BRCA_infor$barcode) ~ "normal"))

tumor_infor <- BRCA_infor %>% filter(tissue == "tumor")
normal_infor <- BRCA_infor %>% filter(tissue == "normal")


tumor_sample <- tcga_counts %>% select(tumor_infor$barcode) %>% as.data.frame()
rownames(tumor_sample) <- rownames(tcga_counts)

normal_sample <- tcga_counts %>% select(normal_infor$barcode) %>% as.data.frame()
rownames(normal_sample) <- rownames(tcga_counts)


##DESeq2
DESeq2_DEG <- function(ecDNA_sample,NonecDNA_sample){
  ##DESeq2
  library("DESeq2")
  DEG_data <- as.data.frame(cbind(ecDNA_sample,NonecDNA_sample)) %>% round(.,digits = 0)
  rownames(DEG_data) <- rownames(tcga_counts)
  group <- as.factor(c(rep("ecDNA",length(ecDNA_sample)), rep("NonecDNA",length(NonecDNA_sample)))) #建立分组
  colGroup <- data.frame(row.names = colnames(DEG_data),
                         group_list = group)
  dds <- DESeqDataSetFromMatrix(countData = DEG_data,
                                colData = colGroup,
                                design = ~ group_list)
  dds <- dds[rowSums(counts(dds)) > 10, ] 
  dds2 <- DESeq(dds)
  res <-  results(dds2, contrast=c("group_list","ecDNA","NonecDNA")) 
  resOrdered <- res[order(res$padj),]
  resOrdered$gene_id <- rownames(resOrdered)
  DE_result <- resOrdered %>% as_data_frame(.) %>% na.omit(.)
  
  return(DE_result)
}
Z_score <- function(data1){
  for (i in 1:length(rownames(data1))) {
    data1 <- as.matrix(data1)
    data1[i,] <- (data1[i,]-mean(data1[i,]))/sd(data1[i,])
  }
  return(data1)
}

DE_Cluster <- DESeq2_DEG(tumor_sample,normal_sample)

fwrite(DE_Cluster,"./model_built/pathway_score_213/time_analysis/bluk_brca_DEG.txt")


###
genelist_four <- KEGG_pathway_gene %>% filter(hsa %in% c("hsa00190","hsa04612","hsa04940","hsa05416"))
genelist_four1 <- genelist_four[,-1]
###

##GSEA
library(GSEABase)
library(clusterProfiler)
HallmarkGeneSet <- read.gmt("./model_built/pathway_score_213/time_analysis/single_gene_list.gmt") 

Gsea_DEG <- DE_Cluster %>% 
  dplyr::mutate(state = (-log10(padj)) * sign(log2FoldChange)) %>% dplyr::arrange(-state) %>%
  dplyr::filter(padj != 0)

geneList <- Gsea_DEG$state 
names(geneList) <- Gsea_DEG$gene_id 
geneList <- sort(geneList, decreasing = T) 
GSEA_result <- GSEA(geneList, TERM2GENE = HallmarkGeneSet, pvalueCutoff = 1,eps = 0) 

library(enrichplot)
gseaplot2(GSEA_result, GSEA_result@result$Description[1:4], title = "", color = "red", base_size = 12,
          rel_heights = c(1.5, 0.5, 1), subplots = 1:3, pvalue_table = T,
          ES_geom = "line")



##GSVA
library(GSVA)
library(GSEABase)
HallmarkGeneSet <- getGmt("./model_built/pathway_score_213/time_analysis/single_gene_list.gmt") 
gsva_result <- gsva(as.matrix(GSE530_sample), HallmarkGeneSet,
                    min.sz=1, max.sz=1000, verbose=TRUE,kcdf="Poisson",parallel.sz=5L)


gsva_result <- as.data.frame(t(gsva_result))

gsva_data1 <- cbind(BRCA_infor,gsva_result)


library(ggpubr)
library(ggprism)
library(ggplot2)
library(cowplot)

gsva_data2 <- melt(gsva_data1)
gsva_data3 <- na.omit(gsva_data2)
gsva_data3$tissue <- factor(gsva_data3$tissue,levels = c("tumor","normal"))


ggplot(data=gsva_data3,aes(x=variable,y=value,fill=factor(tissue)))+
  geom_boxplot()+
  stat_compare_means(aes(label = ..p.signif..))+
  theme_prism()+
  labs(y="GSVA Score",title = "BRCA")+
  theme(axis.title.x = element_blank())




## gene expression

tcga_tpm <- readRDS("~/project/common_data/tcga/tpm_clean_data.rds")
library(NeoEnrichment)
library(dplyr)
library(data.table)
infor <- as.data.frame(colnames(tcga_tpm))
colnames(infor) <- "barcode"
infor$type <- get_cancer_type(infor$barcode)

BRCA_infor <- infor %>% filter(type == "BRCA") 
BRCA_infor <- BRCA_infor %>% dplyr::mutate(tissue = case_when(grepl("*01$",BRCA_infor$barcode) ~ "tumor",
                                                              grepl("*11$",BRCA_infor$barcode) ~ "normal"))

GSE530_sample <- tcga_tpm %>% select(BRCA_infor$barcode) %>% as.data.frame()
rownames(GSE530_sample) <- rownames(tcga_tpm)

data1 <- GSE530_sample[as.character(gene2$Var1),]
data1 <- na.omit(data1)
data2 <- as.data.frame(t(data1))
data3 <- cbind(BRCA_infor, data2)

data4 <- melt(data3)

data4$variable <- factor(data4$variable,levels = c("HLA-A","HLA-B","HLA-C","HLA-E","HLA-F","HLA-G","HLA-DRA",
                                                   "HLA-DRB1","HLA-DRB5","HLA-DQA1","HLA-DQA2",
                                                   "HLA-DQB1","HLA-DOB","HLA-DMA","HLA-DMB","HLA-DOA","HLA-DPA1","HLA-DPB1"))


library(ggpubr)
library(ggprism)
library(ggplot2)
library(cowplot)

data4$tissue <- factor(data4$tissue,levels = c("tumor","normal"))
data4 <- na.omit(data4)

ggplot(data=data4,aes(x=variable,y=value,fill=factor(tissue)))+
  geom_boxplot()+
  stat_compare_means(aes(label = ..p.signif..))+
  theme_prism()+
  theme(axis.text.x = element_text(angle = 30,vjust = 1, hjust = 1) )+
  labs(y="Gene expression",title = "BRCA")+
  theme(axis.title.x = element_blank())



fwrite(gene_list,"./model_built/pathway_score_213/time_analysis/genelist.csv")
