
#' @title data normalized
#' @description Normalize single-cell raw counts matrix.
#' @details Input a data.frame where the rows are the gene names and the columns are the sample names.
#' @param expr_data A single-cell counts expression matrix.
#' @return A normalized single-cell expression matrix.
#' @export
#' @examples



data_normalized <- function(expr_data){

  gene_id <- rownames(expr_data)

  data1 <- expr_data %>% apply(2,function(x){x/sum(x) * 10000}) %>% as.data.frame()
  data2 <- data1 %>% dplyr::mutate_all(funs(log2(.+1)))

  rownames(data2) <- gene_id
  data2 <- round(data2,3)
  return(data2)
}


