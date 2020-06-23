### check zebrafish uORFome for stalling

## get uORFs
libs_path <- "/Volumes/USELESS/META/STALLING/uORFome"
feregg <- get(load(file = c(file.path(libs_path, "leaders_zv10_zf_02_fertilized_egg.csv.uorf.rdata"))))
cell64 <- get(load(file = c(file.path(libs_path, "leaders_zv10_zf_03_64_cells.csv.uorf.rdata"))))
cell512 <- get(load(file = c(file.path(libs_path, "leaders_zv10_zf_04_cells_512.csv.uorf.rdata"))))
rm(rangesOfuORFs)

## get fasta
fasta_cdna <- read.fasta("/Volumes/USELESS/DATA/fasta/danio_rerio/GRCz10/cdna/Danio_rerio.GRCz10.cdna.all.fa")
names(fasta_cdna) <- sapply(names(fasta_cdna), function(x){substr(x,1,18)})

## get ribo coverage (corresponding to libs)

# !!!!!! get GTF for zebrafish corresponding to the files !!!!
# !!!!!! get riboseq too !!!!


