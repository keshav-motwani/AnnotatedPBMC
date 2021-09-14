#' Get 3' PBMC CITE-seq data published in https://www.biorxiv.org/content/10.1101/2020.10.12.335331v1
#'
#' @return
#' @export
get_hao_2020 = function(cache_path) {

  cache_path = file.path(cache_path, "hao_2020")
  dir.create(cache_path)

  return(cache(file.path(cache_path, "hao_2020.rds"),
               function() .process_hao_2020(cache_path)))

}

#' @importFrom BiocFileCache BiocFileCache bfcrpath bfcremove bfcrid removebfc
#' @importFrom Matrix readMM
#' @importFrom SingleCellExperiment altExp altExp<- colData logcounts logcounts<- counts SingleCellExperiment
.process_hao_2020 = function(cache_path) {

  bfc = BiocFileCache(cache_path, ask = FALSE)

  rna_mtx_link = "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5008737&format=file&file=GSM5008737%5FRNA%5F3P%2Dmatrix%2Emtx%2Egz"
  rna_features_link = "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5008737&format=file&file=GSM5008737%5FRNA%5F3P%2Dfeatures%2Etsv%2Egz"
  rna_barcodes_link = "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5008737&format=file&file=GSM5008737%5FRNA%5F3P%2Dbarcodes%2Etsv%2Egz"

  rna_mtx_path = bfcrpath(bfc, rna_mtx_link)
  rna_features_path = bfcrpath(bfc, rna_features_link)
  rna_barcodes_path = bfcrpath(bfc, rna_barcodes_link)

  adt_mtx_link = "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5008738&format=file&file=GSM5008738%5FADT%5F3P%2Dmatrix%2Emtx%2Egz"
  adt_features_link = "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5008738&format=file&file=GSM5008738%5FADT%5F3P%2Dfeatures%2Etsv%2Egz"
  adt_barcodes_link = "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM5008738&format=file&file=GSM5008738%5FADT%5F3P%2Dbarcodes%2Etsv%2Egz"

  adt_mtx_path = bfcrpath(bfc, adt_mtx_link)
  adt_features_path = bfcrpath(bfc, adt_features_link)
  adt_barcodes_path = bfcrpath(bfc, adt_barcodes_link)

  metadata_link = "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE164378&format=file&file=GSE164378%5Fsc%2Emeta%2Edata%5F3P%2Ecsv%2Egz"

  metadata_path = bfcrpath(bfc, metadata_link)

  protein = readMM(adt_mtx_path)
  features = read.table(adt_features_path)
  barcodes = read.table(adt_barcodes_path)

  rownames(protein) = features$V1
  colnames(protein) = barcodes$V1

  gene = readMM(rna_mtx_path)
  features = read.table(rna_features_path)
  barcodes = read.table(rna_barcodes_path)

  rownames(gene) = features$V1
  colnames(gene) = barcodes$V1

  metadata = read.csv(metadata_path)

  sce = SingleCellExperiment(assays = list(counts = gene))
  altExp(sce, "ADT") = SingleCellExperiment(assays = list(counts = protein))

  sce$batch = metadata$Batch
  sce$time = metadata$time
  sce$donor = metadata$donor
  sce$lane = metadata$lane
  sce$cell_type_1 = metadata$celltype.l1
  sce$cell_type_2 = metadata$celltype.l2
  sce$cell_type_3 = metadata$celltype.l3
  sce$phase = metadata$Phase
  sce$dataset = paste0("hao_2020_", sce$donor, sce$time)

  sce = sce[, sce$cell_type_2 != "Doublet"]

  rm(gene, protein)
  gc()

  sce = filter_data(sce)

  logcounts(sce) = normalize_gene(counts(sce))
  logcounts(altExp(sce, "ADT")) = normalize_protein(counts(altExp(sce, "ADT")))

  bfcremove(bfc, bfcrid(bfc))

  return(sce)

}
