library(plyr)

fastaFile <- "upstream.fasta"
bedFile <- "genes.bed"

### Search for motif sequence
motif.search(motif_sequence="taatta", fastaFile=fastaFile, bedFile=bedFile)

motif.search(motif_sequence="taat", fastaFile=fastaFile, bedFile=bedFile)

## Combine the number of forward and reverse motifs in the sequence
motif_sum <- data.frame(motif_search_result$ID)
motif_sum$total <-
  (motif_search_result$fwd_motif_count + motif_search_result$rev_motif_count)
## sort descending
motif_sum <- arrange(motif_sum, desc(total))

### Plot number of motifs for each gene
rows <- as.numeric(nrow(motif_sum))
bar_color <- rainbow(n=rows)
xlimit <- as.numeric(1 + max(motif_sum$total))
par(las=2) # make label text perpendicular to axis
barplot(motif_sum$total, names.arg=motif_sum$motif_search_result.ID,
        horiz=T, col=bar_color, xlab="Number of motifs", ylab="Gene ID",
        xlim=c(0, xlimit), cex.names=0.5)

### Calculate distance between motifs
location_sort <- data.frame(motif_coordinates$ID)
location_sort$start <- motif_coordinates$motifStart
location_sort$end <- motif_coordinates$motifEnd
location_sort <- arrange(location_sort, motif_coordinates.ID, start)

