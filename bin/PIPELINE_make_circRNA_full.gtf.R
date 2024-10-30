#!/usr/bin/env Rscript  

# Function to display the usage of the script  
print_usage <- function() {  
  cat("Usage:\n")  
  cat("Rscript PIPELINE_make_circRNA_full.gtf.R --stout stout.list --outputdir outputdir --reference ReferenceSet.txt\n")  
  cat("\n")  
  cat("Parameters:\n")  
  cat("--stout, -s       : Path to the stout file (required)\n")  
  cat("--outputdir, -d  : Path of the output file (required)\n")  
  cat("--reference, -r   : Path to the reference file (required)\n")  
}  

# Parse command line arguments  
args <- commandArgs(trailingOnly = TRUE)  

# Check the number of parameters  
if (length(args) < 6) {  
  print_usage()  
  stop("Please provide all required parameters.")  
}  

# Define parameters  
stout <- args[which(args == "--stout" | args == "-s") + 1]  
outputdir <- args[which(args == "--outputdir" | args == "-d") + 1]  
reference <- args[which(args == "--reference" | args == "-r") + 1]  

# Load necessary libraries
options(warn = -1)  
suppressPackageStartupMessages({
  library(data.table)
  library(plyr)
  library(dplyr)
  library(rtracklayer)
})

# Start time tracking  
start_time <- Sys.time()  

# Create output directory if it doesn't exist  
output_dir <- dirname(outputdir)  
if (!dir.exists(output_dir)) {  
  dir.create(output_dir, recursive = TRUE)  
} 

# Load reference
ReferenceSet <- fread(reference,data.table = F, sep = "\t", header = TRUE)

# Read input file
stout.list <- fread(stout, data.table = FALSE, sep = "\t", header = FALSE)
# Rename columns with unique names
colnames(stout.list) <- c("Image_ID", "bsj", "chr", "start", "end", "total_exp",
                                      "isoform_number", "isoform_exp", "isoform_length", 
                                      "isoform_state", "strand", "gene_id", "isoform_cirexon")
# Filter for 'Full' isoforms using dplyr
type_full <- stout.list[stout.list$isoform_state == "Full",c('chr', 'start', 'end', 'strand', 'bsj', 'isoform_state', 'isoform_cirexon')]
# Function to generate GTF format for full isoforms
gtf_Full.list = list() 
for(index in 1:nrow(type_full)){
    onerow <- type_full[index,]
    chr <- onerow$chr
    start <- onerow$start
    end <- onerow$end
    strand <- onerow$strand
    isoform_cirexon <- onerow$isoform_cirexon
    exon <- strsplit(strsplit(isoform_cirexon, ",")[[1]], "-")
    exon <- data.frame(t(data.frame(exon)))
    rownames(exon) <- NULL
    colnames(exon) <- c("start","end")
    exon$start <- as.numeric(exon$start)
    exon$end <- as.numeric(exon$end) 
    exonstart <- paste(exon$start, collapse = ",")
    exonend <- paste(exon$end, collapse = ",")
    isoformID <- paste0("chr", chr, "|", exonstart, "|", exonend, "|", strand)
    bsj <- paste0(chr, ":", start, "|", end)
    results <- data.frame(chr = chr, ciri = "ciri", type = "exon", start = exon$start, end = exon$end,
               attr1 = ".", strand = strand, attr2 = ".",bsj = bsj, isoformID = isoformID)
    gtf_Full.list[[index]] <- results           
  }
type_full.gtf <- do.call(rbind, gtf_Full.list)
# final
type_full.gtf$ReferenceSource = mapvalues(type_full.gtf$isoformID,ReferenceSet$isoformID,ReferenceSet$ReferenceSource,warn_missing=F)
type_full.gtf$attr <- paste0('bsj "', type_full.gtf$bsj, '"; ',  
                              'transcript_id "', type_full.gtf$isoformID, '"; ',  
                              'isoform_state "', "Full", '"; ',  
                              'ReferenceSource "', type_full.gtf$ReferenceSource, '"; ')
type_full.gtf = type_full.gtf[,c("chr","ciri","type","start","end","attr1","strand","attr2","attr")]
# Write output to file using fwrite for performance
write.table(type_full.gtf, file = paste0(outputdir,"/circRNA_full.gtf"), sep = "\t", quote = FALSE, col.names = FALSE, append = FALSE, row.names = FALSE)
# End time tracking  
end_time <- Sys.time()  
time_taken <- end_time - start_time
time_taken_minutes <- as.numeric(time_taken, units = "mins")  
# Print the execution time  
print(sprintf("Finished make %s. Time taken: %.2f minutes", paste0(outputdir,"/circRNA_full.gtf"), time_taken_minutes))
