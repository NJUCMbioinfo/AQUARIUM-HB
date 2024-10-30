# script name : PIPELINE_make_reference_set.R  

# Description:  
# This script merges the cirifull data and then generates a reference dataset containing complete circRNA transcript information.  
#  
# Usage:  
# Rscript PIPELINE_make_reference_set.R --input <input_file> --output <output_file>  
# Example: Rscript PIPELINE_make_reference_set.R -i stout.list_path.txt -o test2/ReferenceSet.txt  

# Load necessary libraries  
suppressPackageStartupMessages({  
  library(argparse)    
  library(data.table)  
  library(tidyr)  
  library(plyr)  
  library(dplyr)  
})  

generate_isoform_id <- function(stout_list_full) {
  convert_coordinates <- function(coord_str) {
    coords <- strsplit(coord_str, ",")[[1]]
    start_coords <- sapply(strsplit(coords, "-"), function(x) x[1])
    end_coords <- sapply(strsplit(coords, "-"), function(x) x[2])
    start_str <- paste(start_coords, collapse = ",")
    end_str <- paste(end_coords, collapse = ",")
    return(paste(start_str, "|", end_str, sep = ""))
  }
  stout_list_full$isoformID <- paste0("chr", stout_list_full$chr, "|",
                                       sapply(stout_list_full$isoform_cirexon, convert_coordinates), "|",
                                       stout_list_full$strand)
  return(stout_list_full)
}

# Function to print usage message  
print_usage <- function() {  
  cat("Usage: Rscript PIPELINE_make_reference_set.R --input <input_file> --output <output_file>\n")  
  cat("Example: Rscript PIPELINE_make_reference_set.R -i stout.list_path.txt -o test2/ReferenceSet.txt\n")  
  quit(status = 1)  # Exit after displaying usage  
}  

# Parse arguments from the command line  
parser <- ArgumentParser()  
parser$add_argument("--input", "-i", help = "Input file with list of stout files", required = TRUE)  
parser$add_argument("--output", "-o", help = "Output file for reference set", required = TRUE)  
args <- parser$parse_args()  

# Read the input file (list of file paths and sample IDs)  
filetable <- fread(args$input, data.table = FALSE, header = FALSE)  

# Initialize an empty list to store the merged data  
stout.list_allsample.list <- list()  

# Iterate through each sample and read its respective data file  
for (index in 1:nrow(filetable)) {  
  # Try to read stout list file  
  tryCatch({  
    stout.list <- fread(filetable[index, 1], sep = "\t", data.table = FALSE, header = FALSE)  # Read the stout list file  

    # Check if the data frame is empty  
    if (nrow(stout.list) == 0) {  
      cat("Warning: The file", filetable[index, 1], "is empty and will be skipped.\n")  
      next  # Skip empty files  
    }  

    # Rename the columns to match the expected format  
    names(stout.list) <- c("Image_ID", "bsj", "chr", "start", "end", "total_exp",   
                            "isoform_number", "isoform_exp", "isoform_length",   
                            "isoform_state", "strand", "gene_id", "isoform_cirexon")  

    stout.list <- stout.list[, c("chr", "bsj", "start", "end", "isoform_state", "strand", "isoform_cirexon")]				   

    # Append the current stout.list to the overall data frame  
    stout.list_allsample.list[[index]] <- stout.list   

  }, error = function(e) {  
    # If an error occurs, output the corresponding file name and error message  
    cat("Error: Unable to read the file", filetable[index, 1], "\n")  
    cat("Error message:", conditionMessage(e), "\n")  
  })  
}  

print("Finish importing all stout.list")  

# Combine all samples into one data frame  
stout.list_allsample <- do.call(rbind, stout.list_allsample.list)  
rownames(stout.list_allsample) <- NULL  

print("Filter all full-length circRNAs from stout.list")  
# Filter and remove duplicates  
circRNA_full <- stout.list_allsample[stout.list_allsample$isoform_state == "Full", ]  
circRNA_full <- distinct(circRNA_full)  

circRNA_full <- generate_isoform_id(circRNA_full)  

# Process exon information (Assuming function `parse_transcript` is defined)  
circRNA_full <- circRNA_full %>%  
  rowwise() %>%  
  mutate(  
    parse_results = list(parse_transcript(isoformID)),  
    exon_count = parse_results$exon_count,  
    exon_length = paste(parse_results$exon_length, collapse = ","),  
    exon_total_length = parse_results$exon_total_length  
  ) %>%  
  dplyr::select(-parse_results)   

# Select final required columns  
circRNA_full <- circRNA_full[, c("chr", "bsj", "start", "end", "isoformID", "strand", "exon_count", "exon_length", "exon_total_length")]  

# Create the output directory if it does not exist  
output_dir <- dirname(args$output)  

# Ensure the output directory is not empty and create it if necessary  
if (output_dir != "" && !dir.exists(output_dir)) {  
  dir.create(output_dir, recursive = TRUE)  # Create the directory, including any necessary parent directories  
} else if (output_dir == "") {  
  cat("Error: Output path is not valid.\n")  
  print_usage()  
}  

# Write the full circRNA data to the specified intermediate output path  
temp_output_path <- file.path(output_dir, "circRNAs_from_samples.txt")  
write.table(circRNA_full, file = temp_output_path, sep = "\t", quote = FALSE, col.names = TRUE, append = FALSE, row.names = FALSE)  

# Now we proceed to create the reference set from the generated cirifull output  
cat("Generating reference dataset...\n")  

# Add a column for reference source  
circRNA_full$ReferenceSource <- "Full"  

# Assuming FLcircAS and IsoCirc files are in the current directory or accessible path  
FLcircAS_data <- fread("FLcircAS.final.datatable", data.table = FALSE, header = TRUE)  
IsoCirc_data <- fread("IsoCirc.final.datatable", data.table = FALSE, header = TRUE)  

# Verify if data has been read successfully  
if (nrow(FLcircAS_data) == 0 || nrow(IsoCirc_data) == 0) {  
  cat("Error: FLcircAS or IsoCirc data not found.\n")  
  print_usage()  # Display usage if there's an error  
}  

# Merge all data into a single reference pool  
SupplyReferencePoolraw <- rbind(FLcircAS_data[, names(circRNA_full)],  
                                IsoCirc_data[, names(circRNA_full)],  
                                circRNA_full)  

# Count occurrences of isoformID  
isoform_count <- SupplyReferencePoolraw %>% count(isoformID) %>% arrange(-n)  

# Filter isoform data based on counts  
isoform3 <- SupplyReferencePoolraw[SupplyReferencePoolraw$isoformID %in% isoform_count[isoform_count$n == 3, "isoformID"], ]  
isoform2 <- SupplyReferencePoolraw[SupplyReferencePoolraw$isoformID %in% isoform_count[isoform_count$n == 2, "isoformID"], ]  
isoform1 <- SupplyReferencePoolraw[SupplyReferencePoolraw$isoformID %in% isoform_count[isoform_count$n == 1, "isoformID"], ]  

# Function to summarize isoform data  
summarize_isoform <- function(isoform_data) {  
  isoform_list <- split(isoform_data, isoform_data$isoformID)  
  result_list <- lapply(isoform_list, function(one_isoform) {  
    one_isoform <- one_isoform[order(one_isoform$ReferenceSource, decreasing = FALSE), ]  
    ReferenceSources <- unique(one_isoform$ReferenceSource)  
    ReferenceSource <- paste(ReferenceSources, collapse = ", ")  # Combine reference sources  

    one_isoform$ReferenceSource <- ReferenceSource  
    return(distinct(one_isoform))  # Remove duplicates  
  })  
  return(do.call(rbind, result_list))  
}  

# Summarize isoforms  
isoform3_summary <- summarize_isoform(isoform3)  
isoform2_summary <- summarize_isoform(isoform2)  

# Merge all summarized isoforms into a final reference set  
ReferenceSet <- rbind(isoform1, isoform2_summary, isoform3_summary)  

# Create the output file path for the reference set  
output_path <- args$output  

# Write the final output to the specified file  
write.table(ReferenceSet, file = output_path, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)  

cat("Reference dataset generated and saved to:", output_path, "\n")