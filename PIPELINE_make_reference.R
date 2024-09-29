# Suppress warnings and package startup messages
options(warn = -1)
suppressPackageStartupMessages({
  library(data.table)
  library(tidyr)
  library(plyr)
  library(dplyr)
  library(parallel)
})

isoformID_for_stout.list_Full <- function(stout.list_Full){
  convert_coordinates <- function(coord_str) {
    coords <- strsplit(coord_str, ",")[[1]]
    start_coords <- sapply(strsplit(coords, "-"), function(x) x[1])
    end_coords <- sapply(strsplit(coords, "-"), function(x) x[2])
    start_str <- paste(start_coords, collapse = ",")
    end_str <- paste(end_coords, collapse = ",")
    return(paste(start_str, "|", end_str, sep = ""))
  }
  stout.list_Full$isoformID <- paste0("chr",stout.list_Full$Chr,"|",
                                      sapply(stout.list_Full$isoform_cirexon, convert_coordinates),"|",
                                      stout.list_Full$strain)
  stout.list_Full
}

stout.list <- fread("stout.list_allsample.txt",data.table=F)
stout.list_uniq <- stout.list[,c("Chr","start","end","bsj","isoform_state","strain","isoform_cirexon","isoform_length","gene_id")]
stout.list_uniq <- distinct(stout.list_uniq)
stout.list_full <- stout.list_uniq[stout.list_uniq$isoform_state == "Full",]
stout.list_full <- isoformID_for_stout.list_Full(stout.list_full)
stout.list_full <- stout.list_full[,c("Chr","bsj","strain","isoform_length","gene_id","isoformID")]
names(stout.list_full) <- c("chr","bsj","strand","len","gene_id","isoformID")

dfout_list <- lapply(1:nrow(stout.list_full), function(index) {
  subdata <- stout.list_full[index,]
  subdata <- separate(subdata, bsj, into = c("chr", "start", "end"), remove = FALSE)
  subdata$BSJ_ID <- paste0("chr",subdata$chr,"|",subdata$start,"|",subdata$end,"|",subdata$strand)
  subdata$source <- "Full"
  subdata <- subdata[,c("chr","start","end","strand","BSJ_ID","isoformID","source")]
  return(subdata)
})

isoform_full_uniq.final <- do.call(rbind, dfout_list)
isoform_full_uniq.final$chr <- paste0("chr",isoform_full_uniq.final$chr)


SupplyReferencePoolraw <- rbind(FLcircAS.final.datatable[,names(isoform_full_uniq.final)],
                                IsoCirc.final.datatable[,names(isoform_full_uniq.final)],
                                isoform_full_uniq.final)
names(SupplyReferencePoolraw)
saveRDS(SupplyReferencePoolraw,file=paste0("SupplyReferencePool_raw.rds"))

# 合并同一个转录本不同的 ReferenceType （source）来演： FLcircAS_Blood,FLcircAS_HEK293,FLcircAS_HeLa,FLcircAS_SKOV3,Full,IsoCirc_Blood
args = commandArgs(trailingOnly = TRUE)
chrselect = args[1]
options(stringsAsFactors=F)
options(warn=-1)
library(data.table)
library(tidyr)
library(plyr)
library(dplyr)
dirraw <- "/home/shaoxun/data/circRNA/AQUARIUMfinished/"

SupplyReferencePool <- readRDS(paste0(dirraw,"/FinalResults/SupplyReferencePool_raw.PRJNA722046.rds"))
subdatachr <- SupplyReferencePool[SupplyReferencePool$chr==chrselect,]

print("source for all isoformID")
subdatachr.list <- split(subdatachr,subdatachr$isoformID)
subdatachr.final.list <- list()
for(index in 1:length(subdatachr.list)){
  oneisoform <- subdatachr.list[[index]]
  oneisoform <- oneisoform[order(oneisoform$source,decreasing=FALSE),]
  s1 <- grep("Full",oneisoform$source,value = T)
  s2 <- grep("FLcircAS",oneisoform$source,value = T)
  s3 <- grep("IsoCirc",oneisoform$source,value = T)
  s4 <- grep("Break",oneisoform$source,value = T)
  s5 <- grep("Only",oneisoform$source,value = T)
  all_vars <- list(s1, s2, s3, s4, s5)
  non_empty_vars <- Filter(length, all_vars)
  source <- paste(unlist(non_empty_vars), collapse = ",")
  chr <- unique(oneisoform$chr);start <- unique(oneisoform$start);end <- unique(oneisoform$end)
  strand <- unique(oneisoform$strand);bsj <- paste0(chr,":",start,"|",end)
  bsj <- gsub("chr","",bsj)
  BSJ_ID <- unique(oneisoform$BSJ_ID)
  isoformID <- unique(oneisoform$isoformID)
  exon_start <- strsplit(isoformID, "[|]")[[1]][2]
  exon_end <- strsplit(isoformID, "[|]")[[1]][3]
  oneisoform$exon_start <- exon_start
  oneisoform$exon_end <- exon_end
  len <- sum(as.numeric(unlist(strsplit(exon_end,split = ",")))-as.numeric(unlist(strsplit(exon_start,split = ",")))+1)
  oneisoform <- data.frame(chr=chr,start=start,end=end,strand=strand,len=len,bsj=bsj,BSJ_ID=BSJ_ID,
                           exon_start=exon_start,exon_end=exon_end,isoformID=isoformID,source=source)
  subdatachr.final.list[[index]] <- oneisoform
}
subdatachr.final <- do.call(rbind,subdatachr.final.list)
subdatachr.final <- distinct(subdatachr.final)
outputfile=paste0("/home/shaoxun/data/circRNA/AQUARIUMfinished/FinalResults/SupplyReferencePool_raw.",chrselect,".datatable.PRJNA722046")
write.table(subdatachr.final,file=outputfile,sep="\t",quote=F,col.names=T,append = F,row.names = F)

#data.frame(shell=paste0("nohup Rscript temp.R ",unique(SupplyReferencePoolraw$chr), " &"))
SupplyReferencePool.list <- list()
chrlist <- unique(SupplyReferencePoolraw$chr);chrlist
for(chr in chrlist){
  inputfile=paste0("/home/shaoxun/data/circRNA/AQUARIUMfinished/FinalResults/SupplyReferencePool_raw.",chr,".datatable.PRJNA722046")
  inputdata <- fread(inputfile,data.table = F)
  SupplyReferencePool.list[[chr]] <- inputdata
}
SupplyReferencePool <- do.call(rbind,SupplyReferencePool.list)
names(SupplyReferencePool)[names(SupplyReferencePool)=="isoformID"] <- "transcript_id"
names(SupplyReferencePool)[names(SupplyReferencePool)=="source"] <- "ReferenceType"
write.table(SupplyReferencePool,file=paste0(dirraw,"/FinalResults/SupplyReferencePool.PRJNA722046.txt")  # 这是单数据集的储备池
            ,sep="\t",quote=F,col.names=T,append = F,row.names = F)
