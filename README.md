# AnnotatedPBMC

Provides processed single-cell gene expression and cell type annotations for published peripheral blood mononuclear cell (PBMC) datasets as `Bioconductor SingleCellExperiment` objects. 

The interface is simple: functions are in the form `get_DATASET(cache_path)` where `DATASET` is one of either
  - `10x_pbmc_10k`
  - `10x_pbmc_5k_v3`
  - `10x_sorted`
  - `blish_2020`
  - `ding_2019`
  - `haniffa_2021`
  - `hao_2020`
  - `kotliarov_2020`
  - `su_2020`
  - `tsang_2021`

and `cache_path` is a user-specified file path.

Interested users are encouraged to look at the source code for processing details.

References and original download links are included in the function documentation. Cells are filtered to have at least 200 genes expressed and less than 10% mitochondrial reads, and are log-normalized.
