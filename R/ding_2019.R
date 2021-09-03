#' Get PBMC scRNA-seq data from https://www.biorxiv.org/content/10.1101/632216v2.full.pdf and available https://singlecell.broadinstitute.org/single_cell/study/SCP424/single-cell-comparison-pbmc-data#/
#'
#' @return
#' @export
get_ding_2019 = function(cache_path) {

  cache_path = file.path(cache_path, "ding_2019")
  dir.create(cache_path, recursive = TRUE)

  return(cache(file.path(cache_path, "ding_2019.rds"),
               function() .process_ding_2019(cache_path)))

}

#' @importFrom BiocFileCache BiocFileCache bfcrpath
#' @importFrom Matrix readMM
#' @importFrom SingleCellExperiment altExp altExp<- colData rowData logcounts logcounts<- counts SingleCellExperiment
.process_ding_2019 = function(cache_path) {

  bfc = BiocFileCache(cache_path, ask = FALSE)

  data_link = "https://www.dropbox.com/s/g052oqgjq04c067/counts.umi.txt.gz?dl=1"
  data_path = bfcrpath(bfc, data_link)

  metadata_link = "https://www.dropbox.com/s/lveihz09vbdhrj6/meta.txt?dl=1"
  metadata_path = bfcrpath(bfc, metadata_link)

  genes_link = "https://www.dropbox.com/s/imtmuos4la8crxy/genes.umi.txt?dl=1"
  genes_path = bfcrpath(bfc, genes_link)

  cells_link = "https://www.dropbox.com/s/0fsru2zohj19w9e/cells.umi.new.txt?dl=1"
  cells_path = bfcrpath(bfc, cells_link)

  metadata_names = read.table(metadata_path, sep = "\t", header = FALSE, nrow = 1)
  metadata = read.table(metadata_path, sep = "\t", header = FALSE, skip = 2)
  colnames(metadata) = metadata_names

  genes = read.table(genes_path)
  genes = sapply(strsplit(unlist(genes), "_"), `[`, 2)

  cells = read.table(cells_path)
  cells = unlist(cells)

  data = readMM(data_path)
  rownames(data) = genes
  colnames(data) = cells

  data = data[, intersect(colnames(data), metadata$NAME)]
  names(colnames(data)) = NULL
  names(rownames(data)) = NULL

  metadata = metadata[metadata$NAME %in% intersect(colnames(data), metadata$NAME), ]

  sce = SingleCellExperiment(assays = list(counts = data))

  sce$name = metadata$NAME
  sce$cell_type = metadata$CellType
  sce$experiment = metadata$Experiment
  sce$method = metadata$Method

  sce = sce[, sce$cell_type != "Unassigned"]

  logcounts(sce) = normalize_gene(counts(sce))

  return(sce)

}
