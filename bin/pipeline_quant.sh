#!/bin/bash  
# Usage example:  
# sh pipeline_quant.sh --fastq1 SRR6450118_1.fastq.gz --fastq2 SRR6450118_2.fastq.gz --gtfdir gtfdir --quantdir quantdir  
# 请确保用户将 circumRNA_full.gtf, circRNA_break.gtf, circRNA_only.gtf 文件放在 gtfdir 目录下。  

# 读取配置文件  
source ./configure.cfg  # 确保 fa 和 gtf 的路径在配置文件中定义。  

# 解析输入参数  
while [[ "$#" -gt 0 ]]; do  
    case $1 in  
        --fastq1) fastq1="$2"; shift ;;    # 输入的第一个 FASTQ 文件  
        --fastq2) fastq2="$2"; shift ;;    # 输入的第二个 FASTQ 文件  
        --gtfdir) gtfdir="$2"; shift ;;     # GTF 文件所在目录  
        --quantdir) quantdir="$2"; shift ;; # 量化结果输出目录  
        *) echo "错误：未知参数: $1"; exit 1 ;;   # 处理未知参数  
    esac  
    shift  
done  

# 确保指定目录中有所需的 GTF 文件  
circRNA_full="${gtfdir}/circRNA_full.gtf"  
circRNA_break="${gtfdir}/circRNA_break.gtf"  
circRNA_only="${gtfdir}/circRNA_only.gtf"  

# 检查 quantdir 是否存在，如果不存在则创建它  
if [ ! -d "$quantdir" ]; then  
    mkdir -p "$quantdir"  # 使用 -p 选项以创建父目录  
    echo "已创建目录: $quantdir"  
fi  

# 创建参考线性序列  
gffread "$gtf" -g "$fa" -ME -w "${quantdir}/ref_linear.fa"   

# 合并 circRNA GTF 文件  
cat "$circRNA_full" "$circRNA_break" "$circRNA_only" > "${quantdir}/circRNA_final.gtf"  

# 生成 circRNA 的 raw.fa 文件  
gffread "${quantdir}/circRNA_final.gtf" -g "$fa" -ME -w "${quantdir}/circRNA_raw.fa"  

# 获取 readlen 并运行 R 脚本以添加接头  
readlen=$(zcat "$fastq1" | awk 'NR==2 {print length($0)}')  # 提取 read 长度  
Rscript make.adapt.R "$readlen" "${quantdir}/circRNA_raw.fa" "${quantdir}/circRNA_final.fa"   

# 合并参考和 circRNA 序列  
cat "${quantdir}/ref_linear.fa" "${quantdir}/circRNA_final.fa" > "${quantdir}/final.fa"  
cat "$gtf" "${quantdir}/circRNA_final.gtf" > "${quantdir}/final.gtf"  
    
# 构建 Salmon 索引  
salmon index --kmerLen 31 --transcripts "${quantdir}/final.fa" --index "${quantdir}/index_final" --keepDuplicates  

# 进行量化  
salmon quant \
    --index "${quantdir}/index_final" \
    --libType IU \
    --output "${quantdir}/profile_results" \
    --geneMap "${quantdir}/final.gtf" \
    --mates1 "$fastq1" --mates2 "$fastq2" --threads "$THREAD_COUNT" \
    --seqBias --gcBias --validateMappings  

echo "Finish quant!"