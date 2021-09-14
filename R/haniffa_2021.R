#' Get Healthy PBMC scRNA-seq data from https://www.nature.com/articles/s41591-020-0944-y and available https://www.covid19cellatlas.org/
#'
#' @return
#' @export
get_haniffa_2021 = function(cache_path) {

  cache_path = file.path(cache_path, "haniffa_2021")
  dir.create(cache_path, recursive = TRUE)

  return(cache(file.path(cache_path, "haniffa_2021.rds"),
               function() .process_haniffa_2021(cache_path)))

}

#' @importFrom BiocFileCache BiocFileCache bfcrpath bfcremove bfcrid
#' @importFrom Matrix readMM
#' @importFrom SingleCellExperiment altExp altExp<- colData rowData logcounts logcounts<- counts SingleCellExperiment
#' @importFrom SummarizedExperiment metadata<- assay<- assays<-
#' @importFrom zellkonverter readH5AD
.process_haniffa_2021 = function(cache_path) {

  bfc = BiocFileCache(cache_path, ask = FALSE)

  data_link = "https://covid19.cog.sanger.ac.uk/submissions/release1/haniffa21.processed.h5ad"
  data_path = "/blue/amolstad/keshav.motwani/AnnotatedPBMC/temp_data/haniffa21.processed.h5ad"

  sce = readH5AD(data_path)
  metadata(sce) = list()
  assay(sce, "X") = NULL
  names(assays(sce)) = "counts"
  gc()

  sce = sce[, sce$Status == "Healthy"]
  gc()

  sce$cell_type_1 = sce$initial_clustering
  sce$cell_type_2 = sce$full_clustering
  sce$dataset = paste0("haniffa_2021_", sce$patient_id)

  sce = filter_data(sce)

  logcounts(sce) = normalize_gene(counts(sce))

  bfcremove(bfc, bfcrid(bfc))

  return(sce)

}
