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
