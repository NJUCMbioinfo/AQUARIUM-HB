args = commandArgs(trailingOnly = TRUE)  

# 确保所需参数数量正确  
if (length(args) != 3) {  
    stop("Usage: Rscript script_name.R <adaptlen> <input_fasta_file> <output_fasta_file>")  
}  

adaptlen = as.numeric(args[1])  # 适配器长度，确保转换为数字  
input_fasta_file = args[2]  
output_fasta_file = args[3]  

# 确保输入的适配器长度是正数  
if (is.na(adaptlen) || adaptlen <= 0) {  
    stop("Error: 'adaptlen' must be a positive number.")  
}  

# 设置选项  
options(stringsAsFactors = FALSE)  
options(warn = -1)  

# 加载必要的库  
library(Biostrings)  
library(data.table)  

# 复制输入 FASTA 文件到输出路径  
if (!file.exists(input_fasta_file)) {  
    stop("Error: Input FASTA file does not exist.")  
}  
file.copy(input_fasta_file, output_fasta_file, overwrite = TRUE)  

# 读取输出 FASTA 文件  
output_fasta_seq <- readDNAStringSet(output_fasta_file)  

# 检查是否读取到任何序列  
if (length(output_fasta_seq) == 0) {  
    stop("Error: No sequences found in the input FASTA file.")  
}  

# 为适配器构建新序列  
for (index in 1:length(output_fasta_seq)) {  
    seq = output_fasta_seq[index]  
    
    # 检查序列长度是否有效  
    if (width(seq) == 0) {  
        stop(paste("Error: Sequence", index, "is empty."))  
    }  
    
    # 复制序列以生成适配器  
    seq_for_adaptors <- DNAStringSet(c(rep(unlist(seq), ceiling(adaptlen / width(seq)))))  
    adaptors <- subseq(seq_for_adaptors, start = width(seq_for_adaptors) - adaptlen + 1, end = width(seq_for_adaptors))  

    # 拼接适配器和原始序列  
    new_seqs <- paste0(adaptors, seq)  
    output_fasta_seq[index] <- new_seqs  
}  

# 写入新的 FASTA 文件  
writeXStringSet(output_fasta_seq, output_fasta_file, format = "fasta")