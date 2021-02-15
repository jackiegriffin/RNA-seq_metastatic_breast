# XenofilteR: computational deconvolution of mouse and human reads in tumor xenograft sequence data
# https://github.com/PeeperLab/XenofilteR
# https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-018-2353-5

## install dependencies ----

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("BiocParallel")
BiocManager::install("Rsamtools")
BiocManager::install("GenomicAlignments")
BiocManager::install("GenomicRanges")
BiocManager::install("futile.logger")

library(BiocParallel)
library(Rsamtools)
library(GenomicAlignments)
library(GenomicRanges)
library(futile.logger)
library(IRanges)
library(S4Vectors)

## install XenofilteR package ----

# 1. download tar.gz file from https://github.com/PeeperLab/XenoFilteR/releases
# 2. in terminal, gunzip file 'XenofilteR_1.6.tar.gz'
# 3. move XenofilteR_1.6.tar to desktop R project folder 
# 4. Run the following code to install .tar package:

install.packages("XenofilteR_1.6.tar", repos = NULL, type = "source")
library("XenofilteR")


## generate BAM files for upload ----

# XenofilteR requires only a CIGAR string and NM-tag, standard values present in the BAM-format 

# XenofilteR requires a dataframe or matrix, named 'sample.list',
#      with in the first column the bam file names as mapped to the graft reference. 
#      The second column contains the file names and paths to the bam files as mapped to the host reference.
#      Each row in 'sample.list' represents a single sequence run or sample. 


human_3P <- BamFile('C:/Users/Jgriffin/Desktop/rnaseq_feature/10_S8_human.sorted.bam')
seqinfo(human_3P)
sl_human_3P <- seqlengths(human_3P)
sl_human_3P

quickBamFlagSummary(human_3P)
mouse_3P <- BamFile('C:/Users/Jgriffin/Desktop/rnaseq_feature/3P_mouse_sorted.bam')
seqinfo(mouse_3P)
sl_mouse_3P <- seqlengths(mouse_3P)
sl_mouse_3P
quickBamFlagSummary(mouse_3P)

## run xenofilter ----
# bam files must be sorted based on coordinates; $ samtools sort filename.bam > filename.sorted.bam
sample.list <- read.csv(file = 'sample_list_xeno.csv', stringsAsFactors = FALSE, header = FALSE)



bp.param <- SnowParam(workers = 1, type = "SOCK")
XenofilteR(sample.list, destination.folder = "Xeno_output", bp.param = bp.param)







