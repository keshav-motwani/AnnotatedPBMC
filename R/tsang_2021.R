#' Get Healthy PBMC scRNA-seq data from https://www.nature.com/articles/s41591-020-0944-y and available https://cellxgene.cziscience.com/collections/ed9185e3-5b82-40c7-9824-b2141590c7f0
#'
#' @return
#' @export
get_tsang_2021 = function(cache_path) {

  cache_path = file.path(cache_path, "tsang_2021")
  dir.create(cache_path, recursive = TRUE)

  return(cache(file.path(cache_path, "tsang_2021.rds"),
               function() .process_tsang_2021(cache_path)))

}

#' @importFrom BiocFileCache BiocFileCache bfcrpath bfcremove bfcrid
#' @importFrom Matrix readMM
#' @importFrom SingleCellExperiment altExp altExp<- colData rowData logcounts logcounts<- counts SingleCellExperiment
.process_tsang_2021 = function(cache_path) {

  bfc = BiocFileCache(cache_path, ask = FALSE)

  data_link_1 = "https://corpora-data-prod.s3.amazonaws.com/21d3e683-80a4-4d9b-bc89-ebb2df513dde/local.rds?X-Amz-Algorithm=AWS4-HMAC-SHA256&X-Amz-Credential=ASIATLYQ5N5XXJD6UMVY%2F20210914%2Fus-west-2%2Fs3%2Faws4_request&X-Amz-Date=20210914T132553Z&X-Amz-Expires=604800&X-Amz-SignedHeaders=host&X-Amz-Security-Token=IQoJb3JpZ2luX2VjEOz%2F%2F%2F%2F%2F%2F%2F%2F%2F%2FwEaCXVzLXdlc3QtMiJHMEUCIQCjIQyuaPDjPp0w1HtfPPitzCJl9pwo%2BzAe17dJIi5NOwIgIdnCOncylt%2BE0X5KavnkcKmmeZbxqD1%2BB%2FJrTBRX0%2BUq6wMIRRABGgwyMzE0MjY4NDY1NzUiDCROJItHFPwSciGyIyrIA2PrrguhaMIgKzMs41r7RZ3HyTnipaGneu4Tp766nkq6WMvQNfy0%2FggdQT1Na38xslB5X1ntE6EtNCm2FKRriSy8Arh4DkilqumD%2BNUuWVVxbicaBp7jAczFAh0E3sM0q5wOT99RXtNZFa%2BgoGIN7vE3BfJOs1nYNWTBYaTEIthp4tftCGBcHKsq%2FyATMmls9B5%2B8%2FYELoQ2DUmur2U30EM%2BOyMsoz%2FcroSD9uRGqviIzlfrBvSqL9%2FoJiVBpiSJ%2Bj3%2FyyoUq2RguKtnVFJmFxtgcKebe0ULBMxvXWpU8JvULIGusnZqJH31fWK383R%2BQsNHDT6KiGyclt8vESJIEeVYsvWKN%2B4m%2B%2FaWTgxZibK8l2ryFYVxDE161lq9iqP%2BxB7QtF3vW5FjR0Lmzqh%2BD9coH8KfmAuYDCQJKQ%2Bu%2FiMgzT4rYfVh5GxLFZYc95sfWn5zBcxJI5mkpqDuX4YEvUztnb%2FxwYQtsifvaOCB2U3WonQsL6D%2F0viAGcqnR48PTyv3dwfgwR1Gi%2BFRL28zaGKD%2BB2wVjCBAeM1BmeS1Y5tdDfS4xn2I8T3iPxgCv2Z1W%2FxBG1QcIZUhV8DvZ8k6VaBrKbt%2BTkNIjC4lYKKBjqlASHddeovsajj3GQdnUAlr1mypHcxMdtHGuROWkyMxt2Xvv16xwnd8viCQ4o2GoYYXLAoQrIRbz6%2FcGNqMy2V5ni3u1TVIIEHQi2Z87mMHjstFNp3PZSqVS8Qah9ZPaXMv31lWEJ7JKDHR6E%2Fcr2y%2B16PHyjlflGFdF0SGNBonK1%2FfAgXVKMGx43pMElXk47XMXuCzwhcNicMN9KC7l29TKYNbnY%2BRg%3D%3D&X-Amz-Signature=fa1610e317b26956bb6f70552084ad6117242fce795822234264f65928679bcf"
  data_path_1 = "/blue/amolstad/keshav.motwani/AnnotatedPBMC/temp_data/tsang_adaptive_2021.rds"

  data_link_2 = "https://corpora-data-prod.s3.amazonaws.com/30cd5311-6c09-46c9-94f1-71fe4b91813c/local.rds?X-Amz-Algorithm=AWS4-HMAC-SHA256&X-Amz-Credential=ASIATLYQ5N5XZCONF66P%2F20210914%2Fus-west-2%2Fs3%2Faws4_request&X-Amz-Date=20210914T132625Z&X-Amz-Expires=604800&X-Amz-SignedHeaders=host&X-Amz-Security-Token=IQoJb3JpZ2luX2VjEO3%2F%2F%2F%2F%2F%2F%2F%2F%2F%2FwEaCXVzLXdlc3QtMiJHMEUCIDCQOke0WWS0kx0ltL3EXhaTNk%2FbQOVI7YvfphA8a%2BYyAiEA2iaKA4myW%2FWP2BZHHNfCJyr6t%2BJ9NwYZFc792YEjYMkq6wMIRhABGgwyMzE0MjY4NDY1NzUiDIGiUlwUyGfEQ4Ex3CrIA72UVBLel0QMKqSseNUkoN6A2%2ByYLVa1ipGzTnem3uJNRQZOQv7N40PqLff01ENXI3BT3r25zwgQuaOz8kiN5KhFjnByRhPGruFGBQw61a3jdhIL3AGEA8oPkWml6n%2BAAl2Ctskqb%2BgQCQi9EmXD1V2vAudSkellYrd3AVhAQ4szvLa5leGn2WaeE0pkTVZtWKUAxlM%2B4rBiTwX7ZujyOVXl9Y74cihRENQy4FLnnJqzap48KOsMWRVEzQIEGK4P3tm3GK6Ym1vzzUsiwkTi2YiCZauafryvz53VJzAFghMOs9e8GTgZiL5jynNLjKTvwQMMHwEQ7n%2BOm93Qr1y29PgoRLfddMp9wLN5lp%2BagS%2F%2BwcAmxaHO74f0XL8RXdiWktMHl6M2OMUh0mKPsq6baBMjWm1Mq2kwbAsAsm0BK2NWRRWLqDK4nEkbM5h%2FUIe69nBf2ANv7vDeCtUIPiQF5W9dKOdn9zjObIBWgXwq%2BQKPob3ZTUqiCxVQEpIOjWFyIV6Jjtmsv%2Fucv2HayPLByOjEvWveGhkOhxBXKAELIjBukKi84KxRPrKksWmjRII%2BPeCglermKMrg2RgeyeEz2zEshIWFfWeq3TDyqYKKBjqlARLlO7iSCoKbB8aP05rt4FJP2W2fpb7aq43hRqzuGHmRxcEU5%2BPK6WGajuChluchz0w2IIVBIyz7JNsa8H%2BvpargrgcX5coN3Itucjtrzzezhnq2kqwJterrHWjkK5Gtte0o8RrSdpX2LbQM6GUjI5c0SxGktokAju7LalSuuqvhbj8SOL6jlFqLPffHItxfyQpbrAwwErOrsnbPUHUY2nkpe3KVsQ%3D%3D&X-Amz-Signature=6fcde8cb38f8ba7569be6ea51c86c491cd69f131ae2252ba13ae281d4124846b"
  data_path_2 = "/blue/amolstad/keshav.motwani/AnnotatedPBMC/temp_data/tsang_innate_2021.rds"

  seurat_1 = readRDS(data_path_1)
  seurat_1 = seurat_1[, seurat_1$disease == "Normal"]

  seurat_2 = readRDS(data_path_2)
  seurat_2 = seurat_2[, seurat_2$disease == "Normal"]

  sce = SingleCellExperiment(assays = list(counts = cbind(seurat_1@assays$RNA@counts, seurat_2@assays$RNA@counts)))
  metadata = rbind(seurat_1@meta.data, seurat_2@meta.data)

  rm(seurat_1, seurat_2)
  gc()

  sce$cell_type = metadata$cell_type
  sce$disease = metadata$disease
  sce$donor = metadata$donor
  sce$age = metadata$age
  sce$dataset = paste0("tsang_2021_", sce$donor)

  sce = filter_data(sce)

  logcounts(sce) = normalize_gene(counts(sce))

  bfcremove(bfc, bfcrid(bfc))

  return(sce)

}
