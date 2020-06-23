### search for known stall sites in peaks (which tx?)

### xbp1

xbp1 <- list(yeast = c("YIL101C"),
             fruitfly = c("FBtr0310044", "FBtr0310045"),
             zebrafish = c("ENSDART00000124002", "ENSDART00000124467", "ENSDART00000051649", "ENSDART00000139390"),
             mouse = c("ENSMUST00000063084", "ENSMUST00000145588", "ENSMUST00000149623", "ENSMUST00000149159"),
             human = c("ENST00000216037", "ENST00000403532", "ENST00000484256", "ENST00000405219", "ENST00000344347",
                       "ENST00000482720", "ENST00000611155"))

# termination pausing
# sec61b <- list(zebrafish = c("ENSDART00000059013", "ENSDART00000136980", "ENSDART00000137408", "ENSDART00000147416",
#                              "ENSDART00000146936"),
#                mouse = c("ENSMUST00000065678", "ENSMUST00000137461", "ENSMUST00000125622", "ENSMUST00000136685"),
#                human = c("ENST00000498603", "ENST00000223641", "ENST00000481573"))

peaks_path <- "/Volumes/USELESS/STALLING/peaks_all/NEW_median"
organisms <- c("yeast", "fruitfly", "zebrafish", "mouse", "human")

for (i in 1:length(organisms)) {
  org <- organisms[i]
  libs <- list.files(path = peaks_path, pattern = paste0("^", org)) ## get org-specific libraries
  for (j in 1:length(libs)) {
    load(file = c(file.path(peaks_path, libs[j])))
    print(libs[j])
    print(peaks[peaks$seqnames %in% xbp1[[org]]])
  }
}

# xbp1 comes up in all fruitfly, mouse_3T3, mouse_none - on all tx with z-score above 5.0

# check coverage of xbp1 in gr (is it the selection of well expressed transcripts that excludes it?)

gr_path <- "/Volumes/USELESS/STALLING/gr"

for (i in 1:length(organisms)) {
  org <- organisms[i]
  libs <- list.files(path = gr_path, pattern = paste0("^", org)) ## get org-specific libraries
  for (j in 1:length(libs)) {
    load(file = c(file.path(gr_path, libs[j])))
    gr[seqnames(gr) %in% xbp1[[org]]]
    
    
    print(libs[j])
    
    print(peaks[peaks$seqnames %in% xbp1[[org]]])
  }
}
  
