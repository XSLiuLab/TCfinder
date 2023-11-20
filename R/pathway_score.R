
#' @title pathway score
#' @description Obtain a pathway score matrix for predicting tumor cells.
#' @details Input a sparse matrix, matrix, or data frame where the rows are the gene names and the columns are the sample names. Matrix that can be generated directly using the data_normalized.R function.
#' @param expr_data Single-cell expression matrix after normalization of the original counts matrix.
#' @param normalized If the matrix is not normalized, you need to set normalized = FALSE
#' @param method This parameter is required when normalized = FALSE. If the single-cell sequencing method used is smart-seq2, method = "smart-seq2" is required.
#' For other single-cell sequencing methods, this parameter does not need to be filled in.
#' @param genome This parameter is required when normalized = FALSE. Reference genome, when method = "smart-seq2",
#' this parameter needs to be filled in, you can choose hg19 and hg38
#' @return A matrix containing 213 pathway scores.
#' @export
#' @importFrom dplyr %>% filter select
#' @importFrom Matrix Matrix t
#' @importFrom methods is

pathway_score <- function(expr_data,normalized = TRUE,method = "method",genome = "hg38"){


  if (!methods::is(expr_data, "CsparseMatrix")) {

    expr_data <- Matrix::Matrix(as.matrix(expr_data),sparse = T)

  }


  if (!all(normalized %in% c(TRUE, FALSE))) {
    stop("The normalized parameter is required")
  }


  if(normalized == FALSE){

  normalized_matrix <- TCfinder::data_normalized(expr_data,method = method,genome = genome)

  }

  if(normalized == TRUE){

    normalized_matrix <- expr_data

  }


  KEGG_Gene <- TCfinder::KEGG_Gene
  TCfinder_Pathway <- TCfinder::TCfinder_Pathway

  score_gene <- KEGG_Gene %>% dplyr::filter(hsa %in% TCfinder_Pathway$hsa)


  gene_id <- rownames(normalized_matrix)
  barcode <- colnames(normalized_matrix)
  normalized_matrix <- Matrix::t(normalized_matrix)

  colnames(normalized_matrix) <- gene_id



  all_pathway_score <- NA
  for (i in 1:213) {

    gene <- score_gene %>% dplyr::filter(hsa == names(table(score_gene$hsa))[i])
    pathay_gene <- colnames(normalized_matrix)[which(colnames(normalized_matrix) %in% gene$gene_id==TRUE)]

    selected_data <- normalized_matrix[, pathay_gene]

    path_score <- as.data.frame(apply(selected_data, 1, FUN = function(x){sum(x)/length(x)}))

    colnames(path_score) <- names(table(score_gene$hsa))[i]
    all_pathway_score <- cbind(all_pathway_score,path_score)

  }

  pathway_score <- all_pathway_score[,-1]
  pathway_score <- pathway_score %>% dplyr::select(TCfinder_Pathway$hsa)
  rownames(pathway_score) <- barcode
  return(pathway_score)

}


