#' Get PBMC CITE-seq data published in https://doi.org/10.1038/s41591-020-0769-8
#'
#' @return
#' @export
get_kotliarov_2020 = function(cache_path) {

  cache_path = file.path(cache_path, "kotliarov_2020")
  dir.create(cache_path)

  return(cache(file.path(cache_path, "kotliarov_2020.rds"),
               function() .process_kotliarov_2020(cache_path)))

}

#' @importFrom BiocFileCache BiocFileCache bfcrpath bfcremove bfcrid
#' @importFrom SingleCellExperiment altExp altExp<- colData logcounts logcounts<- counts SingleCellExperiment
.process_kotliarov_2020 = function(cache_path) {

  bfc = BiocFileCache(cache_path, ask = FALSE)

  data_link = "https://nih.figshare.com/ndownloader/files/20706645"
  data_path = bfcrpath(bfc, data_link)

  old_seurat = readRDS(data_path)

  sce = SingleCellExperiment(assays = list(counts = old_seurat@raw.data[, colnames(old_seurat@data)]))
  metadata = old_seurat@meta.data

  cell_type_1_conversion = c(
    C0 = "CD4+ naive T/DNT",
    C1 = "CD4+ memory T",
    C2 = "Monocyte/mDC",
    C3 = "B",
    C4 = "CD8+ memory T",
    C5 = "NK",
    C6 = "CD8+ naive T",
    C7 = "Unconv T",
    C8 = "Non-classical monocyte",
    C9 = "pDC"
  )
  cell_type_2_conversion = c(
    "C0.0.0" = "CD4+ naive T",
    "C0.1.0" = "Double-negative T",
    "C1.0.0" = "CD4+ central and transitional memory T",
    "C1.1.0" = "CD4+ TEMRA and effector memory T",
    "C2.0.0" = "Classical monocytes",
    "C2.0.1" = "IgA+ monocytes",
    "C2.0.2" = "HSC",
    "C2.1.0" = "mDC",
    "C3.0.0" = "Transitional B",
    "C3.0.1" = "Unswitched B",
    "C3.1.0" = "Switched B",
    "C4.0.0" = "CD8+ central and transitional memory T",
    "C4.0.1" = "CD8+ TEMRA and effector memory T",
    "C4.0.2" = "NKT-like",
    "C4.0.3" = "CD8+ CD103+ T",
    "C5.0.0" = "CD16++ NK",
    "C5.0.1" = "CD56hi CD16lo NK",
    "C6.0.0" = "CD8+ naive T",
    "C7.0.0" = "CD161+ double-negative T",
    "C7.0.1" = "Unconventional CD161hi CD8+ T",
    "C8.0.0" = "Non-classical monocytes",
    "C9.0.0" = "pDC"
  )

  sce$cell_type_1 = cell_type_1_conversion[metadata$K1]
  sce$cell_type_2 = cell_type_2_conversion[metadata$K3]

  sce$Barcode = colnames(sce)
  sce$response = metadata$adjmfc.time
  sce$batch = metadata$Batch
  sce$patient = metadata$sampleid

  adt_counts = attributes(old_seurat@assay$CITE)$raw.data

  altExp(sce, "ADT") = SingleCellExperiment(assays = list(counts = adt_counts))

  sce = filter_data(sce)

  logcounts(sce) = normalize_gene(counts(sce))
  logcounts(altExp(sce, "ADT")) = normalize_protein(counts(altExp(sce, "ADT")))

  bfcremove(bfc, bfcrid(bfc))

  return(sce)

}
