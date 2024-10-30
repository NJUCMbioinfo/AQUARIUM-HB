#!/usr/bin/env Rscript  

# Function to display the usage of the script  
print_usage <- function() {  
  cat("Usage example:\n")  
  cat("Rscript PIPELINE_make_circRNA_only.gtf.R --ciri ciri.report --stout stout.list --outputdir circRNA_only.gtf --reference ReferenceSet.txt  \n") 
  cat("Rscript PIPELINE_make_circRNA_only.gtf.R -c ciri.report -s stout.list -d outputdir -r ReferenceSet.txt  \n")   
  cat("\n")  
  cat("Parameters:\n")  
  cat("--ciri, -c            : Path to the circRNA report file (required)\n")  
  cat("--stout, -s           : Path to the stout list file (required)\n")  
  cat("--outputdir, -d      : Path of the output file (required)\n")  
  cat("--reference, -r       : Path to the ReferenceSet set file (required)\n")
  cat("\n")  
  cat("Functionality:\n")  
  cat("This script reads three input files, filters and merges circRNA data,\n")  
  cat("and generates a GTF formatted output file based on exon information.\n")  
}  

# Parse command line arguments  
args <- commandArgs(trailingOnly = TRUE)  

# Check the number of provided parameters  
if (length(args) < 8) {  
  print_usage()  
  stop("Please provide sufficient parameters: --ciri or -c, --stout or -s, --outputfile or -o, --reference or -r ")  
}  

# Define input and output parameters from command line arguments  
ciri <- args[which(args == "--ciri" | args == "-c") + 1]          # Path to the circRNA report file  
stout <- args[which(args == "--stout" | args == "-s") + 1]        # Path to the stout list file   
outputdir <- args[which(args == "--outputdir" | args == "-d") + 1]   # Name of the output file  
reference <- args[which(args == "--reference" | args == "-r") + 1]  # Path to the ReferenceSet set file 

only_supplyfrom_gtf <- function(bsj){
#  print("This function supply only isoform by gtf annotaion")
#  print("inputdata is only bsj")
  chr <- unlist(strsplit(bsj,split = "[:|]"))[1];chr
  start <- as.numeric(unlist(strsplit(bsj,split = "[:|]"))[2]);start
  end <- as.numeric(unlist(strsplit(bsj,split = "[:|]"))[3]);end
  supply_gtf_exon <- gtf_exon[gtf_exon$seqnames == chr & 
                          gtf_exon$start>= as.numeric(start) & gtf_exon$end<= as.numeric(end) ,]   # 使用human gtf
  if(nrow(supply_gtf_exon)>0){
      supply_exon <- supply_gtf_exon[,c("seqnames","start","end")]
      names(supply_exon) <- c("chrom","exonStart","exonEnd")
      supply_exon <- combine_exons(supply_exon)
      supply_exon <- as.data.frame(supply_exon)
      isoformID = paste0("chr",chr,"|",
                             paste(supply_exon$exonStart,collapse = ","),"|",
                             paste(supply_exon$exonEnd,collapse = ","),"|",unique(supply_gtf_exon$strand));
  
      results <- data.frame(chr=chr,start=start,end=end,strand=unique(supply_gtf_exon$strand),bsj=bsj,isoformID=isoformID,
                            isoform_state="onlyoutRef_gtf",ReferenceSource="gtf")
      if (start == supply_exon$exonStart[1] & end == supply_exon$exonEnd[nrow(supply_exon)]) {
          results <- results
      }else{
          results <- data.frame()
      }
  }else{
      results <- data.frame()
  }
}

exontable_to_gtf <- function(exon_df){
#print("This function convert exon to gtf data")
#print("inputfile including 'chr''start''end''strand''bsj''isoformID''isoform_state''ReferenceSource'")   
    gtf.list <- list()
    for(index in 1:nrow(exon_df)){
    onerow <- exon_df[index,]
    chr = onerow$chr
    ciri = "ciri"
    type = "exon" 
    start = onerow$start
    end = onerow$end
    attr1 = "."
    strand = onerow$strand
    attr2 = "."
    attr = paste0('bsj "',onerow$bsj,'"; ','transcript_id "',onerow$isoformID,'"; ',
                  'isoform_state "',onerow$isoform_state,'"; ',
                  'ReferenceSource "',onerow$ReferenceSource,'"; ')
    gtfresults <- data.frame(chr=chr,source="ciri",type="exon",start=start,end=end,
                             attr1=".",strand=strand,attr2=".",attr=attr)
    gtf.list[[index]] <- gtfresults
    }
    gtf <- do.call(rbind,gtf.list)
    gtf                                           
}

# 载入必要的库
options(warn = -1)  
suppressPackageStartupMessages({
  library(data.table)
  library(tidyr)
  library(rtracklayer)
  library(intervals)
  library(plyr)
  library(dplyr)
  library(Rsubread)
  library(GenomicRanges)
  library("SimRVPedigree") # 合并外显子用，需下载并安装
  library("SimRVSequences") # 合并外显子用，需下载并安装
})

# Create output directory if it doesn't exist  
output_dir <- dirname(outputdir)  
if (!dir.exists(output_dir)) {  
  dir.create(output_dir, recursive = TRUE)  
} 

# Start time tracking  
start_time <- Sys.time()  

# 读取参考文件
ReferenceSet <- fread(reference, data.table = FALSE, sep = "\t", header = TRUE)
gtf_exon <- readRDS("/data2/shaoxun/BloodCircleR/dataset/ref/Reference.gtf.rds")
ciri.report <- fread(ciri, data.table = FALSE, sep = "\t", header = TRUE)
stout.list <- fread(stout, data.table = FALSE, sep = "\t", header = FALSE)

# 清洗数据
names(stout.list) <- c("Image_ID", "bsj", "chr", "start", "end", "total_exp",
                        "isoform_number", "isoform_exp", "isoform_length",
                        "isoform_state", "strand", "gene_id", "isoform_cirexon")
stout.list$chr <- gsub("chr", "", stout.list$chr)

# 获取不同类型的 BSJ
type_only <- setdiff(ciri.report$circRNA_ID, stout.list$bsj)
type_onlyinRef <- intersect(ReferenceSet$bsj, type_only)
type_onlyoutRef <- setdiff(type_only, ReferenceSet$bsj)

# 处理来自参考的类型
type_onlyinRef.df <- ReferenceSet[ReferenceSet$bsj %in% type_onlyinRef, ]
type_onlyinRef.df.list <- split(type_onlyinRef.df, type_onlyinRef.df$bsj)
type_onlyinRef.list <- list()

for (index in 1:length(type_onlyinRef.df.list)) {
  onelist <- type_onlyinRef.df.list[[index]]

  if (length(grep("Full|Blood", onelist$ReferenceType)) > 0) {
    select_isoform_for_bsj <- onelist[grep("Full|Blood", onelist$ReferenceSource), ]
    select_isoform_for_bsj <- select_isoform_for_bsj[which.max(select_isoform_for_bsj$exon_total_length), ]
  } else {
    select_isoform_for_bsj <- onelist[which.max(onelist$exon_total_length), ]
  }

  type_onlyinRef.list[[index]] <- select_isoform_for_bsj
}

type_onlyinRef.df <- do.call(rbind, type_onlyinRef.list)
type_onlyinRef.df$isoform_state <- "onlyinRef_ref"
type_onlyinRef.df <- type_onlyinRef.df[, c("chr", "start", "end", "strand", "bsj", "isoformID", "isoform_state", "ReferenceSource")]

# 处理来自 GTF 的类型
if (length(type_onlyoutRef) > 0) {
  type_onlyoutRef.list <- list()
  for (index in 1:length(type_onlyoutRef)) {
    bsj <- type_onlyoutRef[index]
    results <- only_supplyfrom_gtf(bsj)
    type_onlyoutRef.list[[index]] <- results
  }
  
  if (length(type_onlyoutRef.list) > 0) {
    type_onlyoutRef.df <- do.call(rbind, type_onlyoutRef.list)
    type_onlyoutRef.df <- type_onlyoutRef.df[!duplicated(type_onlyoutRef.df), ]
  } else {
    type_onlyoutRef.df <- data.frame()
  }
} else {
  type_onlyoutRef.df <- data.frame()
}

#print("合并所有结果，然后将外显子表格转换为gtf表格")
type_only.df <- rbind(type_onlyinRef.df, type_onlyoutRef.df)

# 外显子的转换为 GTF 格式
type_only_exon.list <- list()
for (index in 1:nrow(type_only.df)) {
  onerow <- type_only.df[index, ]
  exon_start <- as.numeric(unlist(strsplit(unlist(strsplit(onerow$isoformID, split = "\\|"))[[2]], split = ",")))
  exon_end <- as.numeric(unlist(strsplit(unlist(strsplit(onerow$isoformID, split = "\\|"))[[3]], split = ",")))
  
  exon <- data.frame(start = exon_start, end = exon_end)
  exon <- exon[order(exon$start, decreasing = FALSE), ]
  
  multirow <- onerow[rep(1, nrow(exon)), ]
  multirow$start <- exon$start
  multirow$end <- exon$end
  
  type_only_exon.list[[index]] <- multirow
}

type_only_exon <- do.call(rbind, type_only_exon.list)
rownames(type_only_exon) <- NULL
type_only_gtf <- exontable_to_gtf(type_only_exon)
type_only_gtf$chr <- gsub("chr", "", type_only_gtf$chr)

write.table(type_only_gtf, file = paste0(outputdir,"/circRNA_only.gtf"), sep = "\t", quote = FALSE, col.names = FALSE, append = FALSE, row.names = FALSE)

# End time tracking  
end_time <- Sys.time()  
time_taken <- end_time - start_time
time_taken_minutes <- as.numeric(time_taken, units = "mins")  
# Print the execution time  
print(sprintf("Finished make %s. Time taken: %.2f minutes", paste0(outputdir,"/circRNA_only.gtf"), time_taken_minutes))