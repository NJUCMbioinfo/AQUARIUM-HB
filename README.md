<!DOCTYPE html>  
<html lang="zh">  

<body>  

<h1><strong><em>AQUARIUM-HB</em></strong></h1>  
<p>A bioinformatics pipeline to identify, annotate, quantify, and analyze human blood circular RNAs from RNA-seq data.</p>  
<!-- <p>For detailed documentation, please visit the <a href="https://njucmbioinfo.github.io/AQUARIUM-HB/">Project Documentation</a>.</p>  -->

<h2><strong><em>USAGE</em></strong></h2>  

<h3><strong><em>Requirements</em></strong></h3>  
<ul>  
    <li>Linux (CentOS >=7.6 or Ubuntu >=20.04)</li>  
    <li>R (>= 4.1.0)</li>  
    <li>Python (>= 3.8.10)</li>  
    <li>Java (1.8.0_291)</li>  
    <li>Sratoolkit (2.11.3)</li>  
    <li>BWA (0.7.17)</li>  
    <li>Salmon (1.5.0)</li>  
</ul>  

<h3><strong><em>Required R Packages</em></strong></h3>  
<ul>  
    <li>FactoMineR (2.4)</li>  
    <li>clusterProfiler (3.18.1)</li>  
    <li>dplyr (1.0.7)</li>  
    <li>plyr (1.8.6)</li>  
    <li>rtracklayer (1.52.0)</li>  
    <li>msigdbr (7.2.1)</li>  
    <li>tidyr (1.1.3)</li>  
    <li>factoextra (1.0.7)</li>  
</ul> 

<h3><strong><em>Data</em></strong></h3>
<p>You need to download high-throughput RNA-seq data that meets the following criteria:</p>  
<ul>  
    <li>rRNA-depleted</li>  
    <li>Paired-end</li>  
    <li>Blood source</li>  
</ul>  

<p>Following files are required:</p>  
<ul>  
    <li>All circRNAs from FLcircAS database (download link: <a href="https://drive.google.com/file/d/1jjMEzCEEaUaUHrJZLME5O5_8fbsmYU6q/view?usp=drive_link" target="_blank">FLcircAS Download</a>)</li>  
    <li>All circRNAs from IsoCirc database (download link: <a href="https://drive.google.com/file/d/1oLk3MDw4kTZDzO7iA9nSNmDLoKqKxYSx/view?usp=drive_link" target="_blank">IsoCirc Download</a>)</li>  
    <li>Human reference genome (fasta, hg38) (Align index files shound be generated using "bwa index")</li>  
    <li>Human gene annotation (gtf, hg38)</li>  
</ul>

<h3><strong><em>Quick Start<em></strong></h3>  
<ul>  
    <li>  
        <strong>Detect circRNA from RNA-seq data.</strong>  
        <div class="code-container">  
            <code>sh AQUARIUM_HB.sh detect <br>  
            --fastq1 &lt;sample_1.fastq.gz&gt; <br>  
            --fastq2 &lt;sample_2.fastq.gz&gt; <br>  
            --outdir &lt;outdir&gt;</code>  
        </div>  
    </li>  
    <li>  
        <strong>Construct a reference set of human blood full-length circRNAs.</strong>  
        <div class="code-container">  
            <code>sh AQUARIUM_HB.sh reference <br>  
            --stoutlist_path &lt;inputfile&gt; <br>  
            --reference_file &lt;outputfile&gt;</code>  
        </div>  
    </li>  
    <li>  
        <strong>Reconstruct incomplete circRNAs from RNA-seq data.</strong>  
        <div class="code-container">  
            <code>sh AQUARIUM_HB.sh reconstruct <br>  
            --cirireport_file &lt;cirireport&gt; <br>  
            --stoutlist_file &lt;stoutlist&gt; <br>  
            --reference_file &lt;reference&gt; <br>  
            --outputdir &lt;outputdir&gt;</code>  
        </div>  
    </li>  
    <li>  
        <strong>Quantify human blood full-length circRNAs from RNA-seq data.</strong>  
        <div class="code-container">  
            <code>sh AQUARIUM_HB.sh quant <br>  
            --fastq1 &lt;sample_1.fastq.gz&gt; <br>  
            --fastq2 &lt;sample_2.fastq.gz&gt; <br>  
            --gtfdir &lt;gtfdir&gt; <br>  
            --quantdir &lt;quantdir&gt;</code>  
        </div>  
    </li>  
</ul> 

</body>  
</html>

<h3><strong><em>Documentation</em></strong></h3>  
Documentation is available online at https://njucmbioinfo.github.io/AQUARIUM-HB/

<h3><strong><em>Authors</em></strong></h3>  
Authors: Shaoxun Yuan(yuanshaoxun@njucm.edu.cn), Wanjun Gu(wanjungu@njucm.edu.cn)<br><br>  
Maintainer: Shaoxun Yuan

<h3><strong><em>Release Notes</em></strong></h3>  

version 1.0.0：Minor update

<h3><strong><em>Citing AQUARIUM-HB</em></strong></h3>  

Shaoxun Yuan, Xue Bai, Linwei Li, Wanjun Gu. AQUARIUM_HB: a bioinformatics pipeline for human blood circular RNA analysis.
