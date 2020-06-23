### codons in P/A sites

# for each org:
# load logo
# get P/A-site nt
# table

organisms <- c("yeast", "fruitfly", "zebrafish", "mouse", "human")
logo_path <- "/Volumes/USELESS/STALLING/PLOTS/logo_stall/2lib"

for (i in 1:length(organisms)) {
  org <- organisms[i]
  logo <- as.character(read.table(file = c(file.path(logo_path, paste0(org, ".txt"))))$V1)
  p <- sort(table(sapply(logo, function(x){substr(x, 16, 18)})), decreasing = T)
  a <- sort(table(sapply(logo, function(x){substr(x, 19, 21)})), decreasing = T)
  print(org)
  print(p)
  print(a)
}

