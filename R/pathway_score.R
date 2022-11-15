
#' @title pathway score
#' @description Obtain a pathway score matrix for predicting tumor cells.
#' @details Input a data.frame where the rows are the gene names and the columns are the sample names. Matrix that can be generated directly using the data_normalized.R function.
#' @param normalized_matrix Single-cell expression matrix after normalization of the original counts matrix.
#' @return A matrix containing 213 pathway scores.
#' @export
#' @importFrom dplyr %>% filter select


pathway_score <- function(normalized_matrix){

  KEGG_Gene <- TCfinder::KEGG_Gene
  TCfinder_Pathway <- TCfinder::TCfinder_Pathway

  score_gene <- KEGG_Gene %>% filter(hsa %in% TCfinder_Pathway$hsa)


  gene_id <- rownames(normalized_matrix)
  normalized_matrix <- as.data.frame(t(normalized_matrix))
  colnames(normalized_matrix) <- gene_id


  myFun1 <- function(number){

    sum(number)/length(number)

  }


  all_pathway_score <- NA
  for (i in 1:213) {

    gene <- score_gene %>% filter(hsa == names(table(score_gene$hsa))[i])

    data1 <- normalized_matrix %>% select(gene$gene_id[which(gene$gene_id %in% colnames(normalized_matrix)==TRUE)])

    path_score <- as.data.frame(apply(data1, 1, myFun1))
    colnames(path_score) <- names(table(score_gene$hsa))[i]

    all_pathway_score <- cbind(all_pathway_score,path_score)

  }

  pathway_score <- all_pathway_score[,-1]
  return(pathway_score)


}


