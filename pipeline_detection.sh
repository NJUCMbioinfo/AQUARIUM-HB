# Global variables, adjust as needed
dir_detect=/      # Directory to store circRNA identification results
dir_vis=/         # Directory to store circRNA visualization results
THREAD_COUNT=     # Number of threads for script execution
file=             # Sample name

# Global paths, ensure these are correct to avoid errors
scriptdir=/home/user/AQUARIUM_HB  # Directory where the scripts are located
fa=/home/user/Homo_sapiens.GRCh38.dna_sm.chromosomes.fa  # Path to the reference fasta file
gtf=/home/user/Homo_sapiens.GRCh38.94.chr.gtf  # Path to the GTF annotation file

# Generate the align.sam file
bwa mem -T 10 -t $THREAD_COUNT ${fa} ${file}_1.fastq ${file}_2.fastq -o ${dir_detect}/align.sam

# Generate ciri.report and ciri.report.log files
perl ${scriptdir}/CIRI2.pl --in ${dir_detect}/align.sam --out ${dir_detect}/ciri.report --ref_file ${fa} --anno ${gtf} --thread_num $THREAD_COUNT

# Generate files starting with "as"
perl ${scriptdir}/CIRI_AS_v1.2.pl --sam ${dir_detect}/align.sam --ciri ${dir_detect}/ciri.report --out ${dir_detect}/as --ref_file ${fa} --anno ${gtf} --output_all yes

# Generate full_ro1_align.txt and full_ro1.fq files
/usr/bin/java -jar ${scriptdir}/CIRI-full-yzy.jar RO1 -1 ${file}_1.fastq -2 ${file}_2.fastq -o ${dir_detect}/full -t $THREAD_COUNT 

# Generate full_ro1.sam file
bwa mem -T 19 -t $THREAD_COUNT ${fa} ${dir_detect}/full_ro1.fq -o ${dir_detect}/full_ro1.sam

# Generate full_ro2_info.txt and full_ro2.sam files
/usr/bin/java -jar ${scriptdir}/CIRI-full-yzy.jar RO2 -r ${fa} -s ${dir_detect}/full_ro1.sam -l 150 -o ${dir_detect}/full 

# Generate full_merge_circRNA_detail.anno file
/usr/bin/java -jar ${scriptdir}/CIRI-full-yzy.jar Merge -r ${fa} -a ${gtf} -c ${dir_detect}/ciri.report -as ${dir_detect}/as_jav.list -ro ${dir_detect}/full_ro2_info.list -o ${dir_detect}/full

# Visualization
/usr/bin/java -jar  ${scriptdir}/ciri_vis.jar -i ${dir_detect}/full_merge_circRNA_detail.anno -l ${dir_detect}/as_library_length.list -d ${dir_vis} -r ${fa} -min 1
