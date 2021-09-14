#' Get Healthy PBMC scRNA-seq data from https://www.nature.com/articles/s41591-020-0944-y and available https://www.covid19cellatlas.org/
#'
#' @return
#' @export
get_blish_2020 = function(cache_path) {

  cache_path = file.path(cache_path, "blish_2020")
  dir.create(cache_path, recursive = TRUE)

  return(cache(file.path(cache_path, "blish_2020.rds"),
               function() .process_blish_2020(cache_path)))

}

#' @importFrom BiocFileCache BiocFileCache bfcrpath bfcremove bfcrid
#' @importFrom Matrix readMM
#' @importFrom SingleCellExperiment altExp altExp<- colData rowData logcounts logcounts<- counts SingleCellExperiment
.process_blish_2020 = function(cache_path) {

  bfc = BiocFileCache(cache_path, ask = FALSE)

  data_link = "https://hosted-matrices-prod.s3-us-west-2.amazonaws.com/Single_cell_atlas_of_peripheral_immune_response_to_SARS_CoV_2_infection-25/blish_covid.seu.rds"
  data_path = bfcrpath(bfc, data_link)

  seurat = readRDS(data_path)

  sce = SingleCellExperiment(assays = list(counts = seurat@assays$RNA@counts))
  metadata = seurat@meta.data

  rm(seurat)
  gc()

  sce$cell_type_1 = metadata$cell.type.coarse
  sce$cell_type_2 = metadata$cell.type.fine
  sce$status = metadata$Status
  sce$sex = metadata$Sex
  sce$donor = metadata$Donor.full
  sce$dataset = paste0("blish_2020_", sce$donor)

  sce = sce[, sce$status == "Healthy"]

  sce = filter_data(sce)

  logcounts(sce) = normalize_gene(counts(sce))

  bfcremove(bfc, bfcrid(bfc))

  return(sce)

}
