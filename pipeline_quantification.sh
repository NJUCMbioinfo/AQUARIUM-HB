#!/bin/bash  
# -cwd    

# Global variables, adjust as needed
dir_detect=       # Directory to store circRNA identification results
dir_vis=          # Directory to store circRNA visualization results
THREAD_COUNT=     # Number of threads for script execution
file=             # Sample name

# Global paths, adjust as needed
scriptdir=/home/user/AQUARIUM_HB/  # Directory where all scripts and reference files are located
fa=${scriptdir}/Homo_sapiens.GRCh38.dna_sm.chromosomes.fa  # Path to the reference fasta file
gtf=${scriptdir}/Homo_sapiens.GRCh38.94.chr.gtf  # Path to the GTF annotation file

gffread $gtf -g $fa -ME -w ${dir_quant}/ref_linear.fa  

cat ${dir_quant}/circRNA_full.gtf ${dir_quant}/circRNA_break.gtf ${dir_quant}/circRNA_only.gtf > ${dir_quant}/circRNA_final.gtf  
  
gffread ${dir_quant}/circRNA_final.gtf -g $fa -ME -w ${dir_quant}/circRNA_raw.fa  
    
cp /data2/shaoxun/BloodCircleR/AQUARIUM-HB/make.adapt.R .  

Rscript make.adapt.R ${file}  
    
cat ${dir_quant}/ref_linear.fa ${dir_quant}/circRNA_final.fa > ${dir_quant}/final.fa  

cat ${gtf} ${dir_quant}/circRNA_final.gtf > ${dir_quant}/final.gtf  
    
salmon index --kmerLen 31 --transcripts ${dir_quant}/final.fa --index ${dir_quant}/index_final --keepDuplicates  

salmon quant \
        --index ${dir_quant}/index_final \
        --libType IU \
        --output ${dir_quant}/profile_results \
        --geneMap ${dir_quant}/final.gtf \
        --mates1 ${file}_1.fastq.gz --mates2 ${file}_2.fastq.gz --threads ${THREAD_COUNT} \
        --seqBias --gcBias --validateMappings  


