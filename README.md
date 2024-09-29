# AQUARIUM-HB
A pipeline to identify, annotate, quantify and analyze human blood circular RNAs from RNA-seq data

(1) . download Full-Length circRNAs from Nanopore Sequencing Datasets
In this initial step, we aim to gather full-length circular RNAs (circRNAs) for downstream analyses. We utilize resources provided by the FLcircAS and IsoCirc databases:

All circRNA from FLcircAS:
https://drive.google.com/file/d/1jjMEzCEEaUaUHrJZLME5O5_8fbsmYU6q/view?usp=drive_link

All circRNA from IsoCirc:
https://drive.google.com/file/d/1oLk3MDw4kTZDzO7iA9nSNmDLoKqKxYSx/view?usp=drive_link

By integrating these resources, we can obtain a rich dataset of circRNA information that serves as the foundation for subsequent data analysis.

(2) . Detection of circRNAs from RNA-seq Data
In this step, we focus on detecting circRNAs from RNA-seq data, using example datasets downloaded from project ID PRJNA722046, which can be accessed through NCBI SRA.

Example datasets are downloaded from  https://www.ncbi.xyz/Traces/study/?acc=PRJNA722046&o=acc_s%3Aa (project ID : PRJNA722046)

Execute the following command:

$ sh pipeline_detection.sh SRR14240730  
This command automatically runs our detection script, and detailed information about the output can be found in the CIRI Cookbook.

(3) . Construction of Human Blood circRNA Reference
The goal of this step is to construct a comprehensive reference library containing information from full-length circRNAs.

3.1 Integration of Full-Length circRNAs from CIRI Output
The following command integrates full-length circRNAs obtained from CIRI output:

$ Rscript PIPELINE_merge_stout.list.R --datasetname PRJNA722046 --input filepath_stout.list.txt --output stout.list.allsample.txt  
This command effectively merges circRNA data across multiple samples for further analysis.

3.2 Make Human Blood circRNA Reference
Using the script below, we generate the final reference list of circRNAs:

$ Rscript PIPELINE_make_reference.R --input stout.list.allsample.txt --output CircRNA_Reference.txt  
The output file, CircRNA_Reference.txt, will include circRNA information from FLcircAS, IsoCirc, and RNA-seq data.

(4) . Completion of Partial circRNAs
In this step, we aim to fill in missing information for partial circRNAs to ensure data completeness.

Execute the following commands:

$ Rscript make.circRNA_full.gtf.sample.R SRR14240730  
$ Rscript make.circRNA_break.gtf.sample.R SRR14240730  
$ Rscript make.circRNA_only.gtf.sample.R SRR14240730  
These commands generate the relevant data for circRNAs corresponding to each sample, enhancing the accuracy and reliability of our dataset.

(5) . Integration of Blood circRNAs
In this phase, we integrate the extracted circRNA information with the sample data.
Run the command:

$ Rscript make.circRNAgtf_to_isoformID.R filepath_stout.list.txt  
This step ensures that circRNA data from different sources can be effectively combined, setting a solid foundation for subsequent analyses.

(6) . Annotation of Blood circRNAs
Finally, we perform annotation on the integrated circRNA data. This step provides biological relevance to the circRNAs, aiding in understanding their functions and potential as biomarkers.













