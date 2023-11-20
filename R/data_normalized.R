
#' @title data normalized
#' @description Normalize single-cell raw counts matrix.
#' @details Input a data.frame where the rows are the gene names and the columns are the sample names.
#' @param expr_data A single-cell counts expression matrix.
#' @param method If the single-cell sequencing method used is smart-seq2, method = "smart-seq2" is required.
#' For other single-cell sequencing methods, this parameter does not need to be filled in.
#' @param genome Reference genome, when method = "smart-seq2",
#' this parameter needs to be filled in, you can choose hg19 and hg38
#' @return A normalized single-cell expression matrix.
#' @export
#' @importFrom Matrix Matrix Diagonal colSums
#' @importFrom methods is


data_normalized <- function(expr_data,method = "method",genome = "hg38"){


  if (!methods::is(expr_data, "CsparseMatrix")) {

    expr_data <- Matrix::Matrix(as.matrix(expr_data),sparse = T)

  }


  if (method == "method") {


    sparse_data1 <- expr_data %*% Matrix::Diagonal(x = 1 / Matrix::colSums(expr_data)) * 10000

    #nonzero_indices <- which(sparse_data1 != 0, arr.ind = TRUE)
    #sparse_data1[nonzero_indices] <- round(log2(sparse_data1[nonzero_indices] + 1), 3)
    sparse_data1 <- round(log2(sparse_data1 + 1), 3)
    return(sparse_data1)

  }


  if (method == "smart-seq2") {

    if (genome == "hg19") {
      gene_length <- hg19
    }

    if (genome == "hg38") {
      gene_length <- hg38
    }


    colnames(gene_length) <- c("gene_name","Length")


    use_gene_length <- gene_length[gene_length$gene_name %in% rownames(expr_data),]
    gene_names <- use_gene_length$gene_name
    selected_rows <- expr_data[gene_names, ]


    compute_result <- function(x) {
      round((x * 1000 * 1000000) / (use_gene_length[, 2] * sum(x * 1000 / use_gene_length[, 2])), 3)
    }


    result_matrix <- as.data.frame(apply(selected_rows, 2, compute_result))
    colnames(result_matrix) <- colnames(selected_rows)
    rownames(result_matrix) <- rownames(selected_rows)

    sparse_data1 <- Matrix::Matrix(as.matrix(result_matrix),sparse = T)
    #nonzero_indices <- which(sparse_data1 != 0, arr.ind = TRUE)
    #sparse_data1[nonzero_indices] <- round(log2(sparse_data1[nonzero_indices] + 1), 3)
    sparse_data1 <- round(log2(sparse_data1 + 1), 3)

    return(sparse_data1)
  }


  if (!all(method %in% c("method", "smart-seq2"))) {
    stop("Method parameter error ")
  }
}


