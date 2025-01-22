#!/bin/bash  

# Load configuration from config file  
source /home/shaoxun/BloodCircleR/AQUARIUMHB_github/configure.cfg  

# 检查参数个数  
if [ "$#" -lt 2 ]; then  
    echo "Usage: $0 <module> <args>"  
    echo "Modules:"  
    echo "  detect --fastq1 <sample_1.fastq.gz> --fastq2 <sample_2.fastq.gz> --outdir <outdir>"  
    echo "  reference --stoutlist_path <inputfile> --reference_prefix <reference_prefix>"  
    echo "  reconstruct --cirireport_file <cirireport> --stoutlist_file <stoutlist> --reference_file <reference> --outputdir <outputdir>"  
    echo "  quant --fastq1 <sample_1.fastq.gz> --fastq2 <sample_2.fastq.gz> --gtfdir <gtfdir> --quantdir <quantdir>"  
    exit 1  
fi  

# 获取功能模块  
module=$1  
shift # 移除第一个参数  

# 根据功能模块选择不同的操作  
case $module in  
    detect)  
        if [ "$#" -ne 6 ]; then  
            echo "Usage: $0 detect --fastq1 <sample_1.fastq.gz> --fastq2 <sample_2.fastq.gz> --outdir <outdir> "  
            exit 1  
        fi  
        fastq1="$2"  # 第一个输入文件  
        fastq2="$4"  # 第二个输入文件  
		outdir="$6"  # 第二个输入文件  
        sh ${scriptsdir}pipeline_detection.sh --fastq1 "$fastq1" --fastq2 "$fastq2" --outdir "$outdir"  
        ;;   
    
    reference)  
        if [ "$#" -ne 4 ]; then  
            echo "Usage: $0 reference --stoutlist_path <inputfile> --reference_prefix <reference_prefix>"
            exit 1  
        fi  
        stoutlist_path="$2"  
        reference_prefix="$4"     
        Rscript ${scriptsdir}PIPELINE_make_reference_set.R -i "$stoutlist_path" -o "$reference_prefix" 
        ;;  
    
    reconstruct)  
		if [ "$#" -ne 8 ]; then  
            echo "Usage: $0 reconstruct --cirireport_file <cirireport_file> --stoutlist_file <stoutlist_file> --reference_file <reference_file> --outputdir <outputdir>"  
            exit 1  
        fi  
        cirireport_file="$2" 
        stoutlist_file="$4"    
		reference_file="$6"    
		outputdir="$8"    
        Rscript ${scriptsdir}PIPELINE_make_circRNA_full.gtf.R --stout "$stoutlist_file" --reference "$reference_file" --outputdir "$outputdir"  
        Rscript ${scriptsdir}PIPELINE_make_circRNA_break.gtf.R --stout "$stoutlist_file" --reference "$reference_file"  --outputdir "$outputdir"
        Rscript ${scriptsdir}PIPELINE_make_circRNA_only.gtf.R --ciri "$cirireport_file" --stout "$stoutlist_file" --reference "$reference_file"  --outputdir "$outputdir" 
        ;;  
        
    quant)  
        if [ "$#" -ne 8 ]; then  
            echo "Usage: $0 quant --fastq1 <sample_1.fastq.gz> --fastq2 <sample_2.fastq.gz> --gtfdir <gtfdir> --quantdir <quantdir>"  
            exit 1  
        fi  
		fastq1="$2"   # 
        fastq2="$4"   # 
		gtfdir="$6"   # 
		quantdir="$8" #    
        sh ${scriptsdir}pipeline_quant.sh --fastq1 "$fastq1" --fastq2 "$fastq2" --gtfdir "$gtfdir" --quantdir "$quantdir"  
        ;;  
    
    *)  
        echo "Invalid module: $module"  
        echo "Modules:"  
        echo "  detect"  
        echo "  reference"  
        echo "  reconstruct"  
        echo "  quant"  
        exit 1  
        ;;  
esac