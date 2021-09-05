#' Filter data based on percentage mitochondrial and total number of genes expressed
#'
#' @param sce SingleCellExperiment
#'
#' @importFrom SingleCellExperiment counts
#' @importFrom Matrix colSums
#'
#' @return
#' @export
filter_data = function(sce) {

  data = counts(sce)

  total_umi = colSums(data)
  mito = data[grepl("^MT-", rownames(data)),]
  mito_umi = colSums(mito)
  pct_mito = mito_umi / total_umi * 100

  n_genes = colSums(data != 0)

  # at least 200 genes expressed and less than 10% mitochondrial
  keep = which((n_genes > 200) & (pct_mito < 10))

  sce$n_genes = n_genes
  sce$pct_mito = pct_mito

  sce = sce[, keep]

  return(sce)

}
