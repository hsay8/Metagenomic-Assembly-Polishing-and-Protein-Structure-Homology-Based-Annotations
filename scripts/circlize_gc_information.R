#!/usr/bin/env Rscript

# ensure on output directory you have trailing backslash
# make sure the directory is already made as well
# this will create many pdfs
# usage :
# ./circlize_gc_information.R -i genome.fasta -o contig_name/contig_name-gc-information.txt

suppressPackageStartupMessages(library(optparse))

library(optparse)
library(seqinr)

option_list <- list(
    make_option(c("-i","--fasta"), type="character", default=NULL,
            help="fasta file of genome",
            dest="genome_filename"),
    make_option(c("-o","--output"), type="character", default=NULL,
            help="path output directory [default = ./]",
            dest="output_directory")
            )

parser <- OptionParser(usage = "%prog -i genome.fasta -o contig_name/contig_name-gc-information.txt [options]",option_list=option_list)

opt = parse_args(parser)

# keep as default foramt for now.
fasta <- read.fasta(file = opt$genome_filename)  

# get length of sequence
fasta.length <- length(fasta[[1]])   

# get vector of start positions for window size.
starts <- seq(1, (fasta.length), by = 1000)

# get vector of GC content for 1000 base windows
gc.content <- numeric(length(starts))

# iterate through each window, and calculate GC content.
for (i in seq(length(starts))) {

        # last chunk will need special processing

        if (i == length(starts)) {
          # only get sequence to the end of the fasta, otherwise
          # NAs will mess everything up
          chunk <- fasta[[1]][starts[i]:fasta.length]
          chunkGC <- GC(chunk)
          # add GC content for window to gc vector
          gc.content[i] <- chunkGC
        } else {
          # chunks are 0-100,
          chunk <- fasta[[1]][starts[i]:(starts[i]+999)]
          chunkGC <- GC(chunk)
          # add GC content for window to gc vector
          gc.content[i] <- chunkGC
        }
     }


# GC Skew = (G - C)/(G + C)
# cumulative gc skew will be summed, should start at the origin.
# fasta sequence needs to imported as a string in this case.
# re-read the fasta file as a string.
fasta.string <- read.fasta(file = opt$genome_filename, as.string=TRUE)  

# get length of sequence
fasta.length <- nchar(fasta.string[[1]])   

# get vector of start positions for window size.
starts <- seq(1, (fasta.length), by = 1000)

# get vector of GC content for 100 base windows
# create vector the length of the starts
gc.skew <- numeric(length(starts))

# taken from https://stat.ethz.ch/pipermail/r-help/2007-November/147069.html
gcskew <- function(x) {
   if (!is.character(x) || length(x) > 1)
   stop("single string expected")
   tmp <- tolower(s2c(x))
   nC <- sum(tmp == "c")
   nG <- sum(tmp == "g")
   if (nC + nG == 0) return(NA)
   return(100 * (nC - nG)/(nC + nG))
}

# for every starting window
for (i in seq_len(length(starts))) {
  # iterate through the entire genome
  # last iteration needs special treatment.
  if (i == length(starts)) {
    # only go to last base
    gc.skew[i] <- gcskew(substr(fasta.string, starts[i], fasta.length))
  } else {
    gc.skew[i] <- gcskew(substr(fasta.string, starts[i], starts[i] + 1000 - 1))
  }
}

# calculate cumulative gc skew for each window.
# basically the equation is going to be i + i-1.

gc.culm <- numeric(length(starts))

for (i in seq(gc.skew)) {

  # first one needs to be treated specially
  if (i == 1) {
    # only take first column
    gc.culm[i] <- gc.skew[i]
  } else {
    # this works for last one because it looks back, not ahead.
    gc.culm[i] <- sum(gc.skew[1:i])
  }
}

final_gc_info <- data.frame(gc_content = gc.content, gc_skew = gc.skew, gc_culm = gc.culm)

write.table(final_gc_info, file = opt$output_directory, col.names=TRUE, row.names=FALSE, quote = FALSE)
