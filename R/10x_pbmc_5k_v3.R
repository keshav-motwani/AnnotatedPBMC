#' Get PBMC scRNA-seq data published in https://support.10xgenomics.com/single-cell-gene-expression/datasets/3.0.2/5k_pbmc_protein_v3  and annotated by https://www.biorxiv.org/content/10.1101/2021.05.20.445014v1
#'
#' @return
#' @export
get_10x_pbmc_5k_v3 = function(cache_path) {

  cache_path = file.path(cache_path, "10X_pbmc_5k_v3")
  dir.create(cache_path)

  return(cache(file.path(cache_path, "10X_pbmc_5k_v3.rds"),
               function() .process_gottardo_annotated(cache_path, "10X_pbmc_5k_v3", "DCs")))

}

#' @importFrom BiocFileCache BiocFileCache bfcrpath bfcremove bfcrid
#' @importFrom zellkonverter readH5AD
#' @importFrom SingleCellExperiment altExp altExp<- colData logcounts logcounts<- counts SingleCellExperiment
#' @importFrom SummarizedExperiment assays assays<-
.process_gottardo_annotated = function(cache_path, dataset, coarse_only) {

  bfc = BiocFileCache(cache_path, ask = FALSE)

  data_link = paste0("https://fh-pi-gottardo-r-eco-public.s3.amazonaws.com/SingleCellDatasets/", dataset, ".h5ad")
  data_path = bfcrpath(bfc, data_link)

  data = readH5AD(data_path)

  data$cell_type_1 = as.character(data$cell_type_l1)
  data$cell_type_2 = as.character(data$cell_type_l2)
  data$cell_type_l1 = NULL
  data$cell_type_l2 = NULL
  data$cell_type_2 = ifelse(data$cell_type_1 %in% coarse_only, data$cell_type_1, data$cell_type_2)

  data = data[, data$cell_type_2 != "undefined"]

  data$dataset = dataset

  names(assays(data)) = "counts"

  data = filter_data(data)

  logcounts(data) = normalize_gene(counts(data))

  bfcremove(bfc, bfcrid(bfc))

  return(data)

}
