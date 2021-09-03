#' Get PBMC scRNA-seq data published in https://support.10xgenomics.com/single-cell-gene-expression/datasets/3.0.2/5k_pbmc_protein_v3  and annotated by https://www.biorxiv.org/content/10.1101/2021.05.20.445014v1
#'
#' @return
#' @export
get_10x_pbmc_5k_v3 = function(cache_path) {

  cache_path = file.path(cache_path, "10X_pbmc_5k_v3")
  dir.create(cache_path)

  return(cache(file.path(cache_path, "10X_pbmc_5k_v3.qs"),
               function() .process_gottardo_annotated(cache_path, "10X_pbmc_5k_v3", "DCs")))

}
