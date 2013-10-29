###################################################
# motif.search

# function for locating motifs in DNA sequence 
# also searches for motif in reverse orientation
# only reports exact matched of motif sequence and reverse thereof

# ***Rules***
# fasta sequence ID must be exactly the same as gene ID in .bed file
# For genes on "+" strand, .bed file must have "start" coordinate as TSS
# For genes on "-" strand, .bed file must have "end" coordinate as TSS

###################################################

### Load the required R libraries
library("seqinr")
library("Biostrings")
library("plyr")

### Main function for searching for motifs in sequence
motif.search <- function(motif_sequence, fastaFile, bedFile) {
  
  ### set motif as characters
  motif_fwd <- as.character(motif_sequence)
  
  ### get motif length
  motif_len <- as.numeric(nchar(motif_fwd))
  
  ### set reverse motif
  motif_rev <- reverse(motif_fwd)
  
  ## Load the fasta file
  DNAsequences <- read.fasta(file=fastaFile, as.string=T, strip.desc=T)
  
  ## get the length of each sequence
  sequence_length <- data.frame(getName(DNAsequences), (getLength(DNAsequences)))
  colnames(sequence_length) <- c("ID", "length")
  
  ## Load the .bed file and create headers
  GeneBed <- read.table(file=bedFile, header=F, sep="\t")
  ### format GeneBed with meaningful header
  colnames(GeneBed) <- c("chromosome", "geneStart", "geneEnd", "ID", "score", "strand")

  ### Functions to get fwd and rev motif positions of each input sequence
  find.motif_fwd.pos <- function(input){
    words.pos(pattern=motif_fwd, text=input, ignore.case=T)
    }

  find.motif_rev.pos <- function(input){
    words.pos(pattern=motif_rev, text=input, ignore.case=T)
    }

  ## apply find motif functions to each line of input sequence 
  result_motif_fwd <<- lapply(DNAsequences, FUN=find.motif_fwd.pos)
  result_motif_rev <<- lapply(DNAsequences, FUN=find.motif_rev.pos)

  ## Count number of times motif_fwd and motif_rev appears in each sequence
  motif_fwd_count <- t(data.frame(lapply(result_motif_fwd, length)))
  motif_rev_count <- t(data.frame(lapply(result_motif_rev, length)))

  ## Merge gene summary results with ID
  all_motif_count_results <- merge(x=motif_fwd_count, y=motif_rev_count, by="row.names")
  all_motif_count_results$ID <- all_motif_count_results$Row.names
  all_motif_count_results$fwd_motif_count <- all_motif_count_results$V1.x
  all_motif_count_results$rev_motif_count <- all_motif_count_results$V1.y
  all_motif_count_results$Row.names <- NULL
  all_motif_count_results$V1.x <- NULL
  all_motif_count_results$V1.y <- NULL
  motif_search_result <<- all_motif_count_results
  
  ### Get location results into dataframe
  motif_locations_fwd <- ldply(result_motif_fwd, data.frame)
  motif_locations_rev <- ldply(result_motif_rev, data.frame)
  
  ### Manipulate data into a workable format for .bed conversion
  fwd_data <- (t(motif_locations_fwd))
  fwd_sorted <- apply(fwd_data,2,sort, decreasing=T)
  fwd_sorted <- data.frame(t(fwd_sorted))
  fwd_sorted$ID <- fwd_sorted$X1
  fwd_sorted$motifStart <- fwd_sorted$X2
  fwd_sorted$motifStart <- as.character(fwd_sorted$motifStart)
  fwd_sorted$motifStart <- as.numeric(fwd_sorted$motifStart)
  fwd_sorted$X1 <- NULL
  fwd_sorted$X2 <- NULL
  ### calculate end coordinate using motif length
  fwd_sorted$motifEnd <- fwd_sorted$motifStart + motif_len
  ### add motif direction identifer
  fwd_sorted$direction <- "fwd"
  
  rev_data <- (t(motif_locations_rev))
  rev_sorted <- apply(rev_data,2,sort, decreasing=T)
  rev_sorted <- data.frame(t(rev_sorted))
  rev_sorted$ID <- rev_sorted$X1
  rev_sorted$motifStart <- rev_sorted$X2
  rev_sorted$motifStart <- as.character(rev_sorted$motifStart)
  rev_sorted$motifStart <- as.numeric(rev_sorted$motifStart)
  rev_sorted$X1 <- NULL
  rev_sorted$X2 <- NULL
  ### calculate end coordinate using motif length
  rev_sorted$motifEnd <- rev_sorted$motifStart + motif_len
  ### add motif direction identifer
  rev_sorted$direction <- "rev"
  
  ### Concatenate and output results for fwd and rev data
  motif_coordinates <<- rbind(fwd_sorted, rev_sorted)
  
  ### Merge GeneBed data with motif coordinates
  results_merge <- merge(x=motif_coordinates, y=GeneBed, by="ID")
  results_merge <- merge(x=results_merge, y=sequence_length, by="ID")
  result1 <- results_merge
  
  ### calculate genomic coordinates for motifs using transcription start site
  ## note the calculations are different for gene on "+" or "-" strands.
  ## the minus or plus one accounts for the upstream sequence end and TSS being different by 1 base
  result1$intervalStart <-
    ifelse(test=(result1$strand == "+"),
           yes= (((result1$geneStart - 1) - result1$length) + (result1$motifStart)),
           no= (((result1$geneEnd + 1) + (result1$length - result1$motifEnd))))
          
  result1$intervalEnd <- 
    ifelse(test=(result1$strand == "+"),
           yes = (((result1$geneStart - 1) - result1$length) + (result1$motifStart + motif_len)),
           no = (((result1$geneEnd + 1) + (result1$length - result1$motifStart))))
    
  ### get data in format for .bed file output
  bed_output <- data.frame(results_merge$chromosome)
  colnames(bed_output) <- c("chrom")
  bed_output$chromStart <- result1$intervalStart
  bed_output$chromEnd <- result1$intervalEnd
  bed_output$name <- result1$ID
  bed_output$score <- result1$score
  bed_output$strand <- result1$strand
  bed_output$thickStart <- result1$intervalStart
  bed_output$thickEnd <- result1$intervalEnd
  bed_output$itemRGB <- result1$direction
  
  ### Add red and blue colours to .bed output.
  ## Red indicates motif was found in fwd sequence orientation
  ## Blue represents motif was found in rev sequence orientation
  bed_output$itemRGB <- ifelse(test=(result1$direction == "fwd"), yes="255,0,0", no="0,0,255")
    
  ### Create header for .bed file
  track_id <- "track name ="
  track_name <- paste(track_id, motif_sequence, sep=" ")
  description <- "decscription=\"Motif locations\""
  visability <- "visibility=2"
  color <- "itemRgb=\"on\""
  header <- paste(track_name, description, visability, color, sep="\t")
  
  ### create .bed file with header
  write.table(x=header, file=paste(motif_sequence, "_results.bed", sep=""), 
              quote=F, sep="\t", col.names=F, row.names=F)
  
  ### Write results to .bed file with header
  filePattern <- paste(motif_sequence, "_results.bed", sep="")
  write.table(x=bed_output, file=filePattern, append=T, 
              quote=F, sep="\t", col.names=F, row.names=F)
  
}


