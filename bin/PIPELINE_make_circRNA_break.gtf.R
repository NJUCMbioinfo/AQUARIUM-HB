#!/usr/bin/env Rscript  

# Function to display the usage of the script  
print_usage <- function() {  
  cat("Usage example:\n")  
  cat("Rscript PIPELINE_make_circRNA_break.gtf.R -s stout.list -o circRNA_break.gtf -r ReferenceSet.txt \n")  
  cat("\n")  
  cat("Parameters:\n")  
  cat("--stout, -s       : Path to the stout list file (output from ciri_full, required)\n")  
  cat("--outputdir, -d  : Path of the output file (required)\n")  
  cat("--reference, -r   : Path to the reference set file (required)\n")  

  cat("\n")  
}  

# Parse command line arguments  
args <- commandArgs(trailingOnly = TRUE)  

# Check the number of provided parameters  
if (length(args) < 6) {  
  print_usage()  
  stop("Please provide sufficient parameters")  
}  

# Define input and output parameters from command line arguments  
stout <- args[which(args == "--stout" | args == "-s") + 1]  # Path to the stout.list file  
outputdir <- args[which(args == "--outputdir" | args == "-d") + 1] # Name of the output file  
reference <- args[which(args == "--reference" | args == "-r") + 1] # Path to the ReferenceSet file  

# Create output directory if it doesn't exist  
output_dir <- dirname(outputdir)  
if (!dir.exists(outputdir)) {  
  dir.create(output_dir, recursive = TRUE)  
} 

# 加载必要的库
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
  library("SimRVPedigree") # 合并外显子用，下载gz安装
  library("SimRVSequences") # 合并外显子用，下载gz安装
})

#必要的函数
isoformID_to_exon = function(isoformID){
    exon_start <- strsplit(strsplit(isoformID, "[|]")[[1]][2],split = ",")[[1]]
    exon_end <- strsplit(strsplit(isoformID, "[|]")[[1]][3],split = ",")[[1]]
    exon <- paste(paste(exon_start, exon_end, sep = "-"),collapse = ",")
    exon
}
break_supplyfrom_gtf <- function(onerow,isoform_state="breakinRef_gtf"){
#    print("This function supply break isoform by gtf annotaion,isoform_state can be revised")
#    print("inputdata is only onerow for break isoform: 'chr''start''end''strand''bsj''isoform_state''isoform_cirexon'")
    isoform_cirexon_ciri <- unlist(strsplit(onerow$isoform_cirexon, ","))
    exon <- strsplit(unlist(strsplit(isoform_cirexon_ciri, ",")), split = "-")
    exon <- data.frame(t(data.frame(exon)));
    rownames(exon) <- NULL
    colnames(exon) <- c("exonStart", "exonEnd")
    exon$chrom <- onerow$chr
    exon <- exon[, c("chrom", "exonStart", "exonEnd")];
    index0 <- which(isoform_cirexon_ciri == "0-0")
    isoform_cirexon_ciri <- strsplit(isoform_cirexon_ciri, '-')
    before0 <- as.integer(isoform_cirexon_ciri[[index0 - 1]][2]);
    after0 <- as.integer(isoform_cirexon_ciri[[index0 + 1]][1]);
    supply_exon <- gtf_exon[gtf_exon$seqnames == onerow$chr & gtf_exon$start >= before0 & gtf_exon$end <= after0, c("seqnames", "start", "end")]
    supply_exon <- distinct(supply_exon)
    names(supply_exon) <- c("chrom", "exonStart", "exonEnd")
    exon_break_supply <- rbind(exon, supply_exon)
    exon_break_supply$exonStart <- as.numeric(exon_break_supply$exonStart)
    exon_break_supply <- exon_break_supply[!exon_break_supply$exonStart == 0, ]
    exon_break_supply$exonEnd <- as.numeric(exon_break_supply$exonEnd)
    exon_break_supply <- combine_exons(exon_break_supply)
    exon_break_supply$exonStart <- as.numeric(exon_break_supply$exonStart)
    exon_break_supply$exonEnd <- as.numeric(exon_break_supply$exonEnd)
    
    exonstart <- paste(exon_break_supply$exonStart,collapse = ",");
    exonend <- paste(exon_break_supply$exonEnd,collapse = ",");
    chr=onerow$chr;start=onerow$start;end=onerow$end;strand=onerow$strand
    len=sum(as.numeric(exon_break_supply$exonEnd)-as.numeric(exon_break_supply$exonStart)+1)  
    BSJ_ID=paste0("chr",onerow$chr,"|",onerow$start,"|",onerow$end,"|",onerow$strand)
    exon_start=exonstart
    exon_end=exonend
    isoformID <- paste0("chr",chr,"|",exonstart,"|",exonend,"|",strand);
    chr_isoform_cirexon=paste0(onerow$chr,"@",onerow$isoform_cirexon);
    results <- data.frame(chr=chr,start=start,end=end,strand=strand,
                          bsj=onerow$bsj,isoformID=isoformID,
                          isoform_state=isoform_state,ReferenceSource="gtf")
    return(results)
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


# Start time tracking  
start_time <- Sys.time()  

gtf_exon <- readRDS("/data2/shaoxun/BloodCircleR/dataset/ref/Reference.gtf.rds")


# 读取输入文件到数据框
stout.list <- fread(stout, data.table = FALSE, sep = "\t", header = FALSE)
ReferenceSet <- fread(reference, data.table = FALSE, sep = "\t", header = TRUE)

names(stout.list) <- c("Image_ID", "bsj", "chr", "start", "end", "total_exp",
                        "isoform_number", "isoform_exp", "isoform_length",
                        "isoform_state", "strand", "gene_id", "isoform_cirexon")
stout.list <- stout.list[, c("chr", "start", "end", "strand", "bsj", "isoform_state", "isoform_cirexon")]

# 过滤出 'Break' 状态的 isoform 数据
type_break <- stout.list[stout.list$isoform_state == "Break", ]

# 过滤出在参考轨道中的 'Break' isoforms
type_breakinRef <- type_break[type_break$bsj %in% intersect(type_break$bsj, ReferenceSet$bsj), ]

# 获取在参考轨道中的 'Break' isoforms 的数据
ReferenceSet_breakinRef <- ReferenceSet[ReferenceSet$bsj %in% type_breakinRef$bsj, ]
type_breakoutRef_gtf <- type_break[type_break$bsj %in% setdiff(type_break$bsj, ReferenceSet$bsj), ]

# 处理在参考内的 'Break' isoforms
type_breakinRef.list <- list()
for(index in 1:nrow(type_breakinRef)) {
    onerow <- type_breakinRef[index, ]
    isoform_cirexon_ciri <- unlist(strsplit(onerow$isoform_cirexon, ","))
    break_index <- which(isoform_cirexon_ciri == "0-0")
    exon_before_break <- paste(isoform_cirexon_ciri[c(1:(break_index - 1))], collapse = ",")
    exon_after_break <- paste(isoform_cirexon_ciri[c(break_index + 1):length(isoform_cirexon_ciri)], collapse = ",")
    
    ReferenceSet_bsj <- ReferenceSet_breakinRef[ReferenceSet_breakinRef$bsj == onerow$bsj, ]
    ReferenceSet_bsj <- ReferenceSet_bsj[order(ReferenceSet_bsj$exon_total_length, decreasing = TRUE), ]
    select_isoform_for_bsj <- data.frame()
  
    for(isoformindex in 1:nrow(ReferenceSet_bsj)) {
        select_referece <- ReferenceSet_bsj[isoformindex, ]
        select_referece$chr_isoform_cirexon <- paste0(onerow$chr, "@", onerow$isoform_cirexon)
        
        # 获取外显子的信息
        supply_exon <- isoformID_to_exon(select_referece$isoformID)
        
        # 检查外显子是否在断点前后
        is_subset_before <- grepl(paste0("^", exon_before_break), supply_exon, perl = TRUE)  
        is_subset_after <- grepl(paste0(exon_after_break, "$"), supply_exon, perl = TRUE) 
        
        if (is_subset_before & is_subset_after) {
            select_isoform_for_bsj <- rbind(select_isoform_for_bsj, select_referece)
        }
    }

    if(nrow(select_isoform_for_bsj) > 0) {
        select_isoform_for_bsj$isoform_state <- "breakinRef_ref"
        
        # 优先选择 Full/Blood 的转录本，其次是其他组织细胞
        if(length(grep("Full|Blood", select_isoform_for_bsj$ReferenceSource)) > 0) {
            select_isoform_for_bsj <- select_isoform_for_bsj[grep("Full|Blood", select_isoform_for_bsj$ReferenceSource), ] 
            select_isoform_for_bsj <- select_isoform_for_bsj[which.max(select_isoform_for_bsj$exon_total_length), ]
        } else {
            select_isoform_for_bsj <- select_isoform_for_bsj[which.max(select_isoform_for_bsj$exon_total_length), ]
        }
        select_isoform_for_bsj <- select_isoform_for_bsj[, c("chr", "start", "end", "strand", "bsj", "isoformID",
                                                               "isoform_state", "ReferenceSource")]
    } else {
        select_isoform_for_bsj <- break_supplyfrom_gtf(onerow, "breakinRef_gtf")
    }
    type_breakinRef.list[[index]] <- select_isoform_for_bsj
}

type_breakinRef.df <- do.call(rbind, type_breakinRef.list)

# 处理在参考外的 'Break' isoforms，供应来自 GTF
if(nrow(type_breakoutRef_gtf) > 0) {
    type_breakoutRef.list <- list()
    for(index in 1:nrow(type_breakoutRef_gtf)) {
        onerow <- type_breakoutRef_gtf[index, ]
        type_breakoutRef.list[[index]] <- break_supplyfrom_gtf(onerow, "breakoutRef_gtf")
    }
    type_breakoutRef.df <- do.call(rbind, type_breakoutRef.list)
} else {
    type_breakoutRef.df <- data.frame()
}

# 合并结果
type_break.df <- rbind(type_breakinRef.df, type_breakoutRef.df)

# 处理外显子信息，从 isoformID 到外显子
type_break_exon.list <- list()
for(index in 1:nrow(type_break.df)) {
    onerow <- type_break.df[index, ]
    exon_start <- strsplit(onerow$isoformID, split = "[|]")[[1]][2]
    exon_end <- strsplit(onerow$isoformID, split = "[|]")[[1]][3]
    exon_start <- unlist(strsplit(exon_start, split = ","))
    exon_end <- unlist(strsplit(exon_end, split = ","))
    
    exon <- data.frame(start = exon_start, end = exon_end) # 外显子组成数据框
    exon <- exon[order(exon$start, decreasing = FALSE), ]
    
    multirow <- onerow[rep(1, nrow(exon)), ]
    multirow$start <- exon$start
    multirow$end <- exon$end
    multirow$start <- as.numeric(multirow$start)
    multirow$end <- as.numeric(multirow$end)
    
    type_break_exon.list[[index]] <- multirow
}

type_break_exon <- do.call(rbind, type_break_exon.list)
rownames(type_break_exon) <- NULL

# 生成 GTF 格式
type_break_gtf <- exontable_to_gtf(type_break_exon)
type_break_gtf$chr = gsub("chr","",type_break_gtf$chr)

write.table(type_break_gtf, file = paste0(outputdir,"/circRNA_break.gtf"), sep = "\t", quote = FALSE, col.names = FALSE, append = FALSE, row.names = FALSE)

# End time tracking  
end_time <- Sys.time()  
time_taken <- end_time - start_time
time_taken_minutes <- as.numeric(time_taken, units = "mins")  
# Print the execution time  
print(sprintf("Finished make %s. Time taken: %.2f minutes", paste0(outputdir,"/circRNA_break.gtf"), time_taken_minutes))