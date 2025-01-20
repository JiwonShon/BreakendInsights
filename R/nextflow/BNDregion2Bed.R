#!/usr/bin/env Rscript
## #! /usr/bin/R

############# For Rscript execution in Server at cherry ################################
# by Sanghyun Kim, 10 Jun 2024                                              ######
# input: DataPath, aa_barcode & ou                                          ######
# output: bed file for bnd regions                                          ######
########################################################################################

if (!suppressWarnings(require("dplyr"))) install.packages("dplyr")
if (!suppressWarnings(require("GenomicRanges"))) BiocManager::install("GenomicRanges")
if (!suppressWarnings(require("rtracklayer"))) install.packages("rtracklayer")
if (!suppressWarnings(require("plyranges"))) BiocManager::install("plyranges")
if (!suppressWarnings(require("tidyverse"))) install.packages("tidyverse")
if (!suppressWarnings(require("optparse"))) install.packages("optparse")
if (!suppressWarnings(require("biovizBase"))) install.packages("biovizBase")

require(dplyr)
require(plyr)
require(GenomicRanges)
require(rtracklayer)
require(tidyverse)
require(plyranges)
require(optparse)
require(biovizBase)

option_list = list(
  make_option(c("--aa_out_path"), type="character", default=NULL, help="AA output path", metavar="character"),
  make_option(c("--aa_barcode"), type="character", default=NULL, help="aa_barcode or pair_barcode", metavar="character"),
  make_option(c("-o", "--out_dir"), type="character", default=NULL, help="Directory where output will be written", metavar="character"),
  make_option(c("--refbuild"), type="character", default="hg19", help="type your builded references hg19 or hg38", metavar="character"),
  make_option(c("--band_width"), type="integer", default=1000, help="band width size ", metavar="character")
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

DataPath = opt$aa_out_path
aa_barcode = opt$aa_barcode
outputpath  = opt$out_dir
refbuild = opt$refbuild
band_width = opt$band_width


print(paste0("input Path: ", DataPath))
print(paste0("AA barcode: ", aa_barcode))
print(paste0("Bams were aligned with ", refbuild))
print(paste0("Output will be stored in: ", outputpath))
print(paste0("band width: ", band_width))

if (dir.exists(outputpath)) {
  print(paste0("Directory already exists :", outputpath, "Process will be continued"))
} else {
  print(paste0("Directory does not exists :", outputpath, " Process will be continued after creating output path"))
  dir.create(paste0(outputpath), recursive = TRUE)
}




# aa_barcode = "0cea9a43-7170-4273-879f-084cb6e14196__698a6ab3-cc32-4d8b-b438-d175d822590a"
# aa_barcode = "0e48684b-e0b4-4625-9f32-7e0331dfcefe__d3037817-9843-4d60-b64f-af284081d989"
# OutputPath = "~/Dropbox/DataDownload/temp"


# Gappath = "/Users/sanghyun/Dropbox/DataDownload/UCSC_GenomeBrowser/GapLocation_hg19.tsv" # from UCSC table browser
# Gap = data.table::fread(Gappath, sep="\t", header=TRUE)
# Gap_gr = as(Gap[which(Gap$chrom %in% c(paste0("chr", 1:22), "chrX", "chrY"))], "GRanges") # 357 regions

# Trim chromosome lengths
data(ideoCyto, package = "biovizBase")
biovizBase::isIdeogram(ideoCyto$hg19)
seqlengths(ideoCyto$hg19)
print(seqlengths(ideoCyto$hg19))

Length_Hg19 = as.integer(c(249250621, 243199373, 198022430, 191154276, 180915260,
                            171115067, 159138663, 146364022, 141213431, 135534747,
                            135006516, 133851895, 115169878, 107349540, 102531392,
                            90354753,  81195210,  78077248,  59128983,  63025520,  
                            48129895,  51304566, 155270560,  59373566))
names(Length_Hg19) = paste0("chr", c(1:22, "X", "Y"))

##### Function ################

.BDNregionBED <- function(DataPath, aa_barcode, outputpath, band_width) {
  
  half_width = band_width    ## (BND - half_width) 에서 (BND + half_width) 까지 region으로 추출.

  path = dir(path=DataPath, full.names = T, recursive=T, pattern=paste0(aa_barcode, "_AA_results"), include.dirs = T)
  path = paste0(path, "/")
  
  CycleFiles = str_sort(grep("_cycles.txt", dir(path), value = TRUE), numeric=TRUE)   #### list of cycle.txt 
  print(CycleFiles)

  if (!is_empty(CycleFiles)) {
    
    GraphFiles = str_sort(grep("_graph.txt", dir(path), value = TRUE), numeric=TRUE)    #### list of graph.txt
    print(GraphFiles)
    stopifnot(identical(length(CycleFiles),length(GraphFiles))) # [# of cycle] should be same as [# of graph].
    n_amplicon= length(CycleFiles)  
    print(n_amplicon)
    Sample_gr = GRanges(seqnames=NULL,ranges=NULL,strand=NULL)
  
    for (j in 1:n_amplicon) {
      
      print(paste0(aa_barcode,": ", j, "-th amplicon"))
      
      cycle_file_path <- paste0(path, CycleFiles[j]) 
      graph_file_path <- paste0(path, GraphFiles[j])

      ##### Read cycle.txt ##############
      pp<-readLines(cycle_file_path)      
      
      Intervals      <-pp[substr(pp,1,nchar("Interval"))=="Interval"]      ## Get Intervals from cycle.txt
      Intervals       <-ldply(Intervals,function(x) {
        xx <-unlist(strsplit(x,split="\t"))
        data.frame(Interval_raw=x,
                   Interval_index=as.integer(xx[2]),
                   chr= xx[3],
                   start=as.integer(xx[4]),
                   end=as.integer(xx[5]),
                   stringsAsFactors = F)
      })
      print(Intervals)
      Cycle_segs      <-pp[substr(pp,1,nchar("Segment"))=="Segment"]       ## Get decomposed cycles from Cycle.txt
      Cycle_segs      <-ldply(Cycle_segs ,function(x) {
        xx <-unlist(strsplit(x,split="\t"))
        data.frame(segment_raw=x,
                   segment_index = as.integer(xx[2]),
                   segment_chr= xx[3],
                   segment_start=as.integer(xx[4]),
                   segment_end=as.integer(xx[5]),
                   stringsAsFactors = F) })
      print(Cycle_segs)
      ##### Read graph.txt ##############
      qq<-readLines(graph_file_path)  
      
      Graph_segs      <-qq[substr(qq,1,nchar("sequence"))=="sequence"]     ## Get all regions within the amplicon.
      Graph_segs        <-ldply(Graph_segs,function(x) {
        xx <-unlist(strsplit(x,split="\t"))
        chr1 = unlist(strsplit(str_sub(xx[2],1,-2), ":"))[1]
        start = unlist(strsplit(str_sub(xx[2],1,-2), ":"))[2]
        chr2 = unlist(strsplit(str_sub(xx[3],1,-2), ":"))[1]
        end = unlist(strsplit(str_sub(xx[3],1,-2), ":"))[2]
        stopifnot(identical(chr1,chr2))
        data.frame(segment_raw=x,
                   chr = chr1, start = as.integer(start), end = as.integer(end),
                   copyNumber= as.numeric(xx[4]), aveCov=as.numeric(xx[5]), size= as.integer(xx[6]),
                   N_read=as.integer(xx[7]),
                   stringsAsFactors = F)
      })
      
      Cycle_gr<- as(Cycle_segs[,c(1,3:5)], "GRanges")      ##Cycle segments >>> genomic range
      Graph_gr <- as(Graph_segs[,c(2:6)], "GRanges")  ##Graph segments >>> genomic range
      Interval_gr <- as(Intervals[,3:5], "GRanges")   ## Interval in the first line of cycle.txt  >>> genomic range

      rm(Cycle_segs,Graph_segs,Intervals, pp, qq)
      
      #######################################################################
      #### Assign IntervalIndex, ecDNA to the redefined segments ########
      #######################################################################
      
      redefinedGR = GenomicRanges::disjoin(Cycle_gr)
      redefinedGR = plyranges:: bind_ranges(redefinedGR, plyranges::setdiff_ranges(Graph_gr, redefinedGR))
      redefinedGR = sort(redefinedGR, by = ~ seqnames + start + end)
      print(redefinedGR)
      #redefinedGR = redefinedGR %>% mutate(aa_barcode = aa_barcode) 
      mcols(redefinedGR)$aa_barcode <-  aa_barcode
      print(redefinedGR)             # redefinedGR$aa_barcode = aa_barcode               # aa_barcode
      mcols(redefinedGR)$amplicon_idx <- paste0("amplicon", j)
      #redefinedGR = redefinedGR %>% mutate(amplicon_idx = paste0("amplicon", j))  # redefinedGR$amplicon_idx = paste0("amplicon", j)  # amplicond_idx
      mcols(redefinedGR)$IntervalIndex <- subjectHits(findOverlaps(redefinedGR, Interval_gr))
      #redefinedGR = redefinedGR %>% mutate(IntervalIndex = subjectHits(findOverlaps(redefinedGR, Interval_gr)))
                     ##redefinedGR$IntervalIndex = subjectHits(findOverlaps(redefinedGR, Interval_gr))   ### Interval index ###
      
      n_seg <- length(redefinedGR)
      mcols(redefinedGR)$idx <- 1:n_seg
      #redefinedGR = redefinedGR %>% mutate(idx = 1:n_seg)     ## redefinedGR$idx = 1:n_seg               ### index of segment ###
      print(1)
      Sample_gr = plyranges:: bind_ranges(Sample_gr, redefinedGR)
    }
    OverlapCheck = findOverlaps(Sample_gr, Sample_gr)
    stopifnot(identical(queryHits(OverlapCheck), subjectHits(OverlapCheck)))
  } 
  else {Sample_gr = GRanges(seqnames=NULL,ranges=NULL,strand=NULL) }

  seqlevelsStyle(Sample_gr) = "UCSC"
  
  seqlevels(Sample_gr) =str_sort(seqlevels(Sample_gr), numeric=TRUE)
  seqlengths(Sample_gr) = Length_Hg19[names(Length_Hg19) %in% seqlevels(Sample_gr)]
  print(Sample_gr)

  startFlank = GRanges(seqnames=seqnames(Sample_gr), IRanges(start = start(Sample_gr) - half_width, end=start(Sample_gr) + half_width))
  print(startFlank)
  endFlank = GRanges(seqnames=seqnames(Sample_gr), IRanges(start = end(Sample_gr) - half_width, end=end(Sample_gr) + half_width))
  print(endFlank)
  BND_bed = reduce_ranges(plyranges::bind_ranges(startFlank, endFlank))
  print(BND_bed)

  # updated: 2024-07-24 ====
  # check output: 7f130a48-c6c4-40f7-b932-8ab8de632835__078e58c4-ac54-45b3-8a5e-9908be0db216_BNDregion.bed (check: -999)
  seqlevels(BND_bed) = seqlevels(Seqinfo(genome="hg19")[paste0("chr", c(1:22, "X", "Y"))])
  seqinfo(BND_bed) = Seqinfo(genome="hg19")[paste0("chr", c(1:22, "X", "Y"))]
  BND_bed  = trim(BND_bed)
  # BND_bed = unlist(subtract(BND_bed, Gap_gr))
  print(BND_bed)
  # ========================

  BND_bed_df = data.frame(chr=gsub("chr","",seqnames(BND_bed)), start=start(BND_bed), end=end(BND_bed))

  print(BND_bed_df)
  #rtracklayer::export.bed(BND_bed, con = paste0(outputpath, "/", aa_barcode, ".bed" )  # Save as a bed file
  write_tsv(BND_bed_df, paste0(outputpath, "/", aa_barcode, "_BNDregion.bed" ), col_names = FALSE)
  return(BND_bed_df)                                                                       # return as a grange. 필요없으면 comment로.. 
}


######### sever에서는.... ################
# working = "/home/sanghyun/Rworking"
# setwd(working)
# Cohort<- "TCGA"   # Cohort = "PCAWG"    # Cohort = "HMF" 
# DataPath = paste0("/mnt/NAS3/storage/EcDNA_Advanced/data/aaSuite_files/",Cohort,"/aaSuite_somatic_ss/v0.1344.2/GRCh37/minCN4.5/cnsizeMin50000/10X/calls")

# ######## local에서 test 용으로 OK ########
# DataPath = "~/Dropbox/DataDownload/temp"      
# aa_barcode = "0cea9a43-7170-4273-879f-084cb6e14196__698a6ab3-cc32-4d8b-b438-d175d822590a"
# aa_barcode = "0e48684b-e0b4-4625-9f32-7e0331dfcefe__d3037817-9843-4d60-b64f-af284081d989"
# OutputPath = "~/Dropbox/DataDownload/temp"


out =.BDNregionBED(DataPath, aa_barcode, outputpath, band_width)
out
