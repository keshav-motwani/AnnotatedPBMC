#' Get sorted subsets of PBMC data published in https://www.nature.com/articles/ncomms14049
#'
#' @return
#' @export
get_10x_sorted = function(cache_path) {

  cache_path = file.path(cache_path, "10x_sorted")
  dir.create(cache_path, recursive = TRUE)

  return(cache(file.path(cache_path, "10x_sorted.rds"),
               function() .process_10x_sorted(cache_path)))

}

#' @importFrom BiocFileCache BiocFileCache bfcrpath bfcremove bfcrid
#' @importFrom DropletUtils read10xCounts
#' @importFrom SingleCellExperiment altExp altExp<- colData rowData logcounts logcounts<- counts SingleCellExperiment
.process_10x_sorted = function(cache_path) {

  labels = c(
    "cd14_monocytes" = "CD14+ Monocytes",
    "b_cells" = "CD19+ B cells",
    "cd34" = "CD34+ Cells",
    "regulatory_t" = "CD4+/CD25+ Regulatory T Cells",
    "naive_t" = "CD4+/CD45RA+/CD25- Naive T cells",
    "memory_t" = "CD4+/CD45RO+ Memory T Cells",
    "cd56_nk" = "CD56+ Natural Killer Cells",
    "naive_cytotoxic" = "CD8+/CD45RA+ Naive Cytotoxic T Cells"
  )

  bfc = BiocFileCache(cache_path, ask = FALSE)

  X_list = vector("list", length(labels))
  names(X_list) = names(labels)
  Y_list = vector("list", length(labels))
  names(Y_list) = names(labels)

  for (cell_type in names(labels)) {

    link = paste0("https://cf.10xgenomics.com/samples/cell-exp/1.1.0/", cell_type, "/", cell_type, "_filtered_gene_bc_matrices.tar.gz")
    path = bfcrpath(bfc, link)

    temp_path = file.path(cache_path, "temp")
    untar(path, exdir = temp_path)

    sce = read10xCounts(paste0(temp_path, "/filtered_matrices_mex/hg19"))
    unlink(temp_path, recursive = TRUE)

    rownames(sce) = rowData(sce)$Symbol
    colnames(sce) = sce$Barcode
    X = counts(sce)
    Y = rep(labels[cell_type], ncol(X))

    X_list[[cell_type]] = X
    Y_list[[cell_type]] = Y

  }

  X = do.call(cbind, X_list)
  Y = unlist(Y_list)

  sce = SingleCellExperiment(assays = list(counts = X))
  sce$cell_type = Y

  sce = filter_data(sce)

  logcounts(sce) = normalize_gene(counts(sce))

  bfcremove(bfc, bfcrid(bfc))

  return(sce)

}
