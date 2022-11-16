
#' @title Cell types prediction.
#' @description Classify tumor cells from normal cells.
#' @details Input the pathway score matrix calculated by the pathway_score function.
#' @param path_score The pathway score matrix calculated by the pathway_score function.
#' @return A data.frame containing cell types and predicted values.
#' @export
#' @importFrom reticulate source_python



predict_cell <- function(path_score){

  barcode <- rownames(path_score)
  Path <- fs::path_package("extdata",package = "TCfinder")
  reticulate::source_python(paste0(Path,"/predict_py.py"))

  predict <- predict_py(path_score,Path)
  predict_result <- as.data.frame(predict)
  result <- predict_result %>% mutate(cell_type = case_when(V1 > 0.5 ~ "normal",
                                                            V1 <= 0.5 ~ "tumor"))
  colnames(result) <- c("value","cell_type")
  result$barcode <- barcode
  return(result)
}
