
#' @title data normalized
#' @description Normalize single-cell raw counts matrix.
#' @details Input a data.frame where the rows are the gene names and the columns are the sample names.
#' @param expr_data A single-cell counts expression matrix.
#' @param method Single cell sequencing type, "10X" or "smart-seq2"
#' @param genome Reference genome, when method == "smart-seq2",
#' this parameter needs to be filled in, you can choose hg19 and hg38
#' @return A normalized single-cell expression matrix.
#' @export
#' @importFrom dplyr %>% filter select mutate_all funs


data_normalized <- function(expr_data,method,genome = "hg38"){


  if (method == "10X") {

    gene_id <- rownames(expr_data)

    data1 <- expr_data %>% apply(2,function(x){x/sum(x) * 10000}) %>% as.data.frame()
    data2 <- data1 %>% dplyr::mutate_all(funs(log2(.+1)))

    rownames(data2) <- gene_id
    data2 <- round(data2,3)
    return(data2)
  }

  if (method == "smart-seq2") {

    if (genome == "hg19") {
      gene_length <- hg19
    }

    if (genome == "hg38") {
      gene_length <- hg38
    }


    colnames(gene_length) <- c("gene_name","Length")
    expr_data$gene_name <- rownames(expr_data)
    use_data <- dplyr::inner_join(gene_length,expr_data) %>% as.data.frame()


    result_value <- use_data
    for (i in 3:ncol(use_data)) {
      result <- round((use_data[,i]*1000*1000000)/(use_data[,2]*sum((use_data[,i]*1000/use_data[,2]))),3)

      result_value[,i] <- result
    }


    data1 <- result_value %>% dplyr::select(-Length,-gene_name) %>% as.data.frame()
    data2 <- data1 %>% dplyr::mutate_all(funs(log2(.+1)))
    rownames(data2) <- use_data$gene_name
    return(data2)
  }



  if (!all(method %in% c("10X", "smart-seq2"))) {
    stop("Please choose a sequencing method 10X or smart-seq2")
  }

}


