#!/bin/bash  

# Load configuration from config file  
source ./configure.cfg  

# Check if the necessary arguments are provided  
if [ "$#" -ne 4 ]; then  
    echo "Usage: $0 -1 <fastq1> -2 <fastq2>"  
    echo "Example:"  
    echo "  $0 -1 sample_1.fastq -2 sample_2.fastq"  
    echo "  $0 -1 sample_1.fastq.gz -2 sample_2.fastq.gz"  
    exit 1  
fi  

# Parse input fastq files from command line arguments  
while getopts "1:2:" opt; do  
    case ${opt} in  
        1 )  
            fastq1=$OPTARG  
            ;;  
        2 )  
            fastq2=$OPTARG  
            ;;  
        \? )  
            echo "Usage: $0 -1 <fastq1> -2 <fastq2>"  
            echo "Example:"  
            echo "  $0 -1 sample_1.fastq -2 sample_2.fastq"  
            echo "  $0 -1 sample_1.fastq.gz -2 sample_2.fastq.gz"  
            exit 1  
            ;;  
    esac  
done  

# Extract sample name from fastq1 (assuming the file is named like sample_1.fastq or sample_1.fastq.gz)  
file=$(basename "$fastq1" _1.fastq.gz)  
file=$(basename "$file" _1.fastq)  

# Define directories for results  
dir_detect="./${file}/full"     # Directory to store circRNA identification results  
dir_vis="./${file}/vis"         # Directory to store circRNA visualization results  

# Create necessary directories if they don't exist  
mkdir -p "$dir_detect" "$dir_vis"  

# Generate the align.sam file  
${path_bwa} mem -T 10 -t $THREAD_COUNT ${fa} "$fastq1" "$fastq2" -o "${dir_detect}/align.sam"  

# Generate ciri.report and ciri.report.log files  
${path_perl} ${scriptdir}/CIRI2.pl --in "${dir_detect}/align.sam" --out "${dir_detect}/ciri.report" --ref_file ${fa} --anno ${gtf} --thread_num $THREAD_COUNT  

# Generate files starting with "as"  
${path_perl} ${scriptdir}/CIRI-AS.pl --sam "${dir_detect}/align.sam" --ciri "${dir_detect}/ciri.report" --out "${dir_detect}/as" --ref_file ${fa} --anno ${gtf} --output_all yes  

# Generate full_ro1_align.txt and full_ro1.fq files  
${path_java} -jar ${scriptdir}/CIRI-full.jar RO1 -1 "$fastq1" -2 "$fastq2" -o "${dir_detect}/full" -t $THREAD_COUNT   

# Generate full_ro1.sam file  
${path_bwa} mem -T 19 -t $THREAD_COUNT ${fa} "${dir_detect}/full_ro1.fq" -o "${dir_detect}/full_ro1.sam"  

# Generate full_ro2_info.txt and full_ro2.sam files  
${path_java} -jar ${scriptdir}/CIRI-full.jar RO2 -r ${fa} -s "${dir_detect}/full_ro1.sam" -l 150 -o "${dir_detect}/full"   

# Generate full_merge_circRNA_detail.anno file  
${path_java} -jar ${scriptdir}/CIRI-full.jar Merge -r ${fa} -a ${gtf} -c "${dir_detect}/ciri.report" -as "${dir_detect}/as_jav.list" -ro "${dir_detect}/full_ro2_info.list" -o "${dir_detect}/full"  

# Visualization  
${path_java} -jar  ${scriptdir}/CIRI-vis.jar -i "${dir_detect}/full_merge_circRNA_detail.anno" -l "${dir_detect}/as_library_length.list" -d "${dir_vis}" -r ${fa} -min 1