#!/bin/bash  

source ./configure.cfg  

# 检查所需参数的数量是否正确  
if [ "$#" -ne 6 ]; then  
    echo "Usage: $0 --fastq1 <fastq1> --fastq2 <fastq2> --outdir <outdir>"  
    echo "Example:"  
    echo "  $0 --fastq1 sample_1.fastq.gz --fastq2 sample_2.fastq.gz --outdir sample"  
    exit 1  
fi  

# 解析输入参数  
while [[ "$#" -gt 0 ]]; do  
    case $1 in  
        --fastq1) fastq1="$2"; shift ;;    # 输入的第一个 FASTQ 文件  
        --fastq2) fastq2="$2"; shift ;;    # 输入的第二个 FASTQ 文件  
        --outdir) outdir="$2"; shift ;;     # 输出目录  
        *) echo "未知参数: $1"; exit 1 ;;    # 处理未知参数  
    esac  
    shift  
done  

# 定义结果目录  
dir_detect="${outdir}/full"  # 存储 circRNA 识别结果的目录  
dir_vis="${outdir}/vis"        # 存储 circRNA 可视化结果的目录  

# 创建必要的目录（如果不存在）  
mkdir -p "$dir_detect" "$dir_vis"  

# 生成 align.sam 文件  
${path_bwa} mem -T 10 -t $THREAD_COUNT ${fa} "$fastq1" "$fastq2" -o "${dir_detect}/align.sam"  

# 生成 ciri.report 和 ciri.report.log 文件  
${path_perl} ${scriptsdir}CIRI2.pl --in "${dir_detect}/align.sam" --out "${dir_detect}/ciri.report" --ref_file ${fa} --anno ${gtf} --thread_num $THREAD_COUNT  

# 生成以 "as" 开头的文件  
${path_perl} ${scriptsdir}CIRI-AS.pl --sam "${dir_detect}/align.sam" --ciri "${dir_detect}/ciri.report" --out "${dir_detect}/as" --ref_file ${fa} --anno ${gtf} --output_all yes  

# 生成 full_ro1_align.txt 和 full_ro1.fq 文件  
${path_java} -jar ${scriptsdir}CIRI-full.jar RO1 -1 "$fastq1" -2 "$fastq2" -o "${dir_detect}/full" -t $THREAD_COUNT   

# 生成 full_ro1.sam 文件  
${path_bwa} mem -T 19 -t $THREAD_COUNT ${fa} "${dir_detect}/full_ro1.fq" -o "${dir_detect}/full_ro1.sam"  

# 生成 full_ro2_info.txt 和 full_ro2.sam 文件  
${path_java} -jar ${scriptsdir}CIRI-full.jar RO2 -r ${fa} -s "${dir_detect}/full_ro1.sam" -l 150 -o "${dir_detect}/full"   

# 生成 full_merge_circRNA_detail.anno 文件  
${path_java} -jar ${scriptsdir}CIRI-full.jar Merge -r ${fa} -a ${gtf} -c "${dir_detect}/ciri.report" -as "${dir_detect}/as_jav.list" -ro "${dir_detect}/full_ro2_info.list" -o "${dir_detect}/full"  

# 可视化  
${path_java} -jar ${scriptsdir}CIRI-vis.jar -i "${dir_detect}/full_merge_circRNA_detail.anno" -l "${dir_detect}/as_library_length.list" -d "${dir_vis}" -r ${fa} -min 1  

echo "Finished!"