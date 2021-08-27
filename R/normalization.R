#' Normalize protein expression by scaling each feature by geometric mean
#'
#' @param data matrix with protein expression data
#'
#' @importFrom Matrix rowMeans Diagonal
#'
#' @return
#' @export
normalize_protein = function(data) {

  feature_scales = exp(rowMeans(log1p(data)))

  normalized = log1p(Diagonal(x = 1/feature_scales) %*% data)

  return(normalized)

}

#' Normalize gene expression to log counts per 10,000
#'
#' @param data matrix with protein expression data
#'
#' @importFrom Matrix colSums Diagonal
#'
#' @return
#' @export
normalize_gene = function(data) {

  size_factors = colSums(data)

  data = log1p(data %*% Diagonal(x = 10000/size_factors))

  return(data)

}
