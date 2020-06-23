### plot read/offset
library(data.table)
library(ggplot2)

read_count_path <- "/Volumes/USELESS/STALLING/PLOTS/logo_stall/read_count"


## yeast
y1 <- read.table(file = c(file.path(read_count_path, "yeast_brar.txt")))
y2 <- read.table(file = c(file.path(read_count_path, "yeast_chx.txt")))
y3 <- read.table(file = c(file.path(read_count_path, "yeast_none.txt")))

offsets_y1 <- data.table(read_length = c(21:34), offsets = rep(12, 14))
offsets_y2 <- data.table(read_length = c(24:30), offsets = c(8:14))
offsets_y3 <- data.table(read_length = c(25:32), offsets = c(9:14, 14, 14))

reads_offsets <- list(r = list(y1, y2, y3), o = list(offsets_y1, offsets_y2, offsets_y3))

## fruit fly
f1 <- read.table(file = c(file.path(read_count_path, "fruitfly_dunn.txt")))
f2 <- read.table(file = c(file.path(read_count_path, "fruitfly_luo_DMSO.txt")))
f3 <- read.table(file = c(file.path(read_count_path, "fruitfly_SRR.txt")))

offsets_f1 <- data.table(read_length = c(26:40), offsets = c(rep(12,5),13,rep(14,8),15))
offsets_f2 <- data.table(read_length = c(22:34), offsets = c(rep(8,3),rep(11,2),rep(12,2),rep(13,3),rep(14,3)))
offsets_f3 <- data.table(read_length = c(21:35), offsets = c(rep(5,3),6,7,9,c(9:13),rep(13,4)))

reads_offsets <- list(r = list(f1, f2, f3), o = list(offsets_f1, offsets_f2, offsets_f3))

## zebrafish
z1 <- read.table(file = c(file.path(read_count_path, "zebrafish_bazzini.txt")))
z2 <- read.table(file = c(file.path(read_count_path, "zebrafish_beaudoin.txt")))
z3 <- read.table(file = c(file.path(read_count_path, "zebrafish_subtelny2h.txt")))
z4 <- read.table(file = c(file.path(read_count_path, "zebrafish_dome.txt")))
z5 <- read.table(file = c(file.path(read_count_path, "zebrafish_subtelny4h.txt")))
z6 <- read.table(file = c(file.path(read_count_path, "zebrafish_shield.txt")))
z7 <- read.table(file = c(file.path(read_count_path, "zebrafish_subtelny6h.txt")))

offsets_z1 <- data.table(read_length = c(20:31), offsets = c(3:10,rep(12,3),13))
offsets_z2 <- data.table(read_length = c(25:31), offsets = c(9:12,12:14))
offsets_z3 <- data.table(read_length = c(26,27,32), offsets = rep(12,3))
offsets_z4 <- data.table(read_length = c(25:31), offsets = c(8,9,11,11,12,12,13))
offsets_z5 <- data.table(read_length = c(30), offsets = c(12))
offsets_z6 <- data.table(read_length = c(25:32), offsets = c(8,8,11,11,12,12,13,13))
offsets_z7 <- data.table(read_length = c(30), offsets = c(12))

reads_offsets <- list(r = list(z1, z2, z3, z4, z5, z6, z7),
                      o = list(offsets_z1, offsets_z2, offsets_z3, offsets_z4, offsets_z5, offsets_z6, offsets_z7))

## mouse
m1 <- read.table(file = c(file.path(read_count_path, "mouse_chx.txt")))
m2 <- read.table(file = c(file.path(read_count_path, "mouse_none.txt")))
m3 <- read.table(file = c(file.path(read_count_path, "mouse_3T3.txt")))

offsets_m1 <- data.table(read_length = c(21:22,27:35), offsets = c(rep(3,2),rep(12,4),13,rep(14,4)))
offsets_m2 <- data.table(read_length = c(25:35,40), offsets = c(7:14,rep(14,4)))
offsets_m3 <- data.table(read_length = c(24:35,40), offsets = c(rep(11,3),12,rep(13,8),14))

reads_offsets <- list(r = list(m1, m2, m3), o = list(offsets_m1, offsets_m2, offsets_m3))

## human
h1 <- read.table(file = c(file.path(read_count_path, "human_stern_chx.txt")))
h2 <- read.table(file = c(file.path(read_count_path, "human_stern_none.txt")))
h3 <- read.table(file = c(file.path(read_count_path, "human_stumpf.txt")))
h4 <- read.table(file = c(file.path(read_count_path, "human_subtelny.txt")))

offsets_h1 <- data.table(read_length = c(22:36,40), offsets = c(rep(5,3),rep(8,3),11,12,12,13,13,rep(14,5)))
offsets_h2 <- data.table(read_length = c(21:35,40), offsets = c(3:14,rep(14,4)))
offsets_h3 <- data.table(read_length = c(20:35), offsets = c(3:6,6:9,rep(12,3),13,13,rep(14,3)))
offsets_h4 <- data.table(read_length = c(25:35), offsets = rep(13,11))

reads_offsets <- list(r = list(h1, h2, h3, h4), o = list(offsets_h1, offsets_h2, offsets_h3, offsets_h4))

##########

end5 <- data.table(offsets = c(3:16), reads = rep(0,14))
end3 <- data.table(offsets = c(12:25), reads = rep(0,14))

for (n in range(1:length(reads_offsets[[1]]))) {
  libX <- reads_offsets$r[[n]]
  offsetsX <- reads_offsets$o[[n]]
  for (i in 1:nrow(offsetsX)) {
    end5$reads[which(end5$offsets == offsetsX[i,]$offsets)] <-
      end5$reads[which(end5$offsets == offsetsX[i,]$offsets)] + libX[libX$V2 == offsetsX[i,]$read_length,]$V1
    end3$reads[which(end3$offsets == (offsetsX[i,]$read_length - offsetsX[i,]$offsets))] <- 
      end3$reads[which(end3$offsets == (offsetsX[i,]$read_length - offsetsX[i,]$offsets))] +
      libX[libX$V2 == offsetsX[i,]$read_length,]$V1
  }
}

#f <- data.table(offset = c(8:14), reads = c(158337,	521197,	920331,	3800338, 39761964,	13828310,	8413981))
ggplot(end5, aes(x = offsets, y = reads)) + geom_bar(stat = "identity") + scale_x_reverse() + theme_void()
ggplot(end3, aes(x = offsets, y = reads)) + geom_bar(stat = "identity") + theme_void()


