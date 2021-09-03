#' Get PBMC CITE-seq data published in https://www.sciencedirect.com/science/article/pii/S0092867420314446 and annotated by https://www.biorxiv.org/content/10.1101/2021.05.20.445014v1
#'
#' @return
#' @export
get_su_2020 = function(cache_path) {

  cache_path = file.path(cache_path, "su_2020")
  dir.create(cache_path)

  return(cache(file.path(cache_path, "su_2020.qs"),
               function() .process_gottardo_annotated(cache_path, "su_2020", c())))

}

#' @importFrom BiocFileCache BiocFileCache bfcrpath
#' @importFrom zellkonverter readH5AD
#' @importFrom SingleCellExperiment altExp altExp<- colData logcounts logcounts<- counts SingleCellExperiment
#' @importFrom SummarizedExperiment assays assays<-
.process_gottardo_annotated = function(cache_path, dataset, coarse_only) {

  bfc = BiocFileCache(cache_path, ask = FALSE)

  data_link = paste0("https://fh-pi-gottardo-r-eco-public.s3.amazonaws.com/SingleCellDatasets/", dataset, ".h5ad")
  data_path = bfcrpath(bfc, data_link)

  data = readH5AD(data_path)

  data$cell_type_l1 = as.character(data$cell_type_l1)
  data$cell_type_l2 = as.character(data$cell_type_l2)
  data$cell_type_l2 = ifelse(data$cell_type_l1 %in% coarse_only, data$cell_type_l1, data$cell_type_l2)

  data = data[, data$cell_type_l2 != "undefined"]

  data$cell_type = data$cell_type_l2

  names(assays(data)) = "counts"

  logcounts(data) = normalize_gene(counts(data))

  return(data)

}
