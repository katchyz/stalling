# convert gr to dt
library(data.table)

granges_path <- "/Volumes/USELESS/DATA/Ribo-Seq/fantom_human/granges"
datatable_path <- "/Volumes/USELESS/DATA/Ribo-Seq/fantom_human/datatable"
libs <- list.files(path = granges_path)

for (i in 1:length(libs)) {
  ## load file
  load(file = c(file.path(granges_path, libs[i]))) # gr
  dt <- as.data.frame(gr)
  setDT(dt)
  save(dt, file = c(file.path(datatable_path, libs[i])))
}

