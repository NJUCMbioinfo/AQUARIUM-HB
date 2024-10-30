<!DOCTYPE html>  
<html lang="zh">  

<body>  

<h1><strong><em>AQUARIUM-HB</em></strong></h1>  
<p>A bioinformatics pipeline to identify, annotate, quantify, and analyze human blood circular RNAs from RNA-seq data.</p>  
<p>For detailed documentation, please visit the <a href="https://njucmbioinfo.github.io/AQUARIUM-HB/">Project Documentation</a>.</p>  

<h2><strong><em>USAGE</em></strong></h2>  

<h3><strong><em>Requirements</em></strong></h3>  
<ul>  
    <li>Linux</li>  
    <li>R</li>  
    <li>Python</li>  
    <li>Java</li>  
    <li>Sratoolkit</li>  
    <li>BWA</li>  
    <li>Salmon</li>  
</ul>  

<h3><strong><em>Required R Packages</em></strong></h3>  
<ul>  
    <li>FactoMineR</li>  
    <li>clusterProfiler</li>  
    <li>dplyr</li>  
    <li>plyr</li>  
    <li>rtracklayer</li>  
    <li>msigdbr</li>  
    <li>tidyr</li>  
    <li>factoextra</li>  
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
    <li>Human reference genome (fasta, hg38) (Index files shound be generated using "bwa index *.fasta")</li>  
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
