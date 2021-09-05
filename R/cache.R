cache = function(file_path, fn) {

  dir.create(dirname(file_path), recursive = TRUE)

  if (grepl(".csv", file_path)) {
    write = readr::write_csv
    read = readr::read_csv
  } else {
    write = function(object, path)
      saveRDS(object,  path)
    read = function(path)
      readRDS(path)
  }

  if (!file.exists(file_path)) {
    result = fn()
    write(result, file_path)
  } else {
    result = read(file_path)
  }

  return(result)

}
