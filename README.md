<h1><strong><em>AQUARIUM-HB</em></strong></h1>

A bioinformatics pipeline to identify, annotate, quantify, and analyze human blood circular RNAs from RNA-seq data.

For detailed documentation, please visit the [Project Documentation](https://njucmbioinfo.github.io/AQUARIUM-HB/).

<h1><strong><em>USAGE</em></strong></h1>

<h2><strong><em>Requirements</em></strong></h2>  
<ul>  
    <li>Linux</li>  
    <li>R</li>  
    <li>Python</li>  
    <li>Java</li>  
    <li>Sratoolkit</li>  
    <li>BWA</li>  
    <li>Salmon</li>  
</ul>  

<h3><strong>Required R Packages</strong></h3>  
<p>The following R packages are required:</p>  
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

<h2><strong><em>Data</em></strong></h2>
<p>You need to download high-throughput RNA-seq data that meets the following criteria:</p>
<ul>
    <li>rRNA-depleted</li>
    <li>Paired-end</li>
    <li>Blood source</li>
</ul>

<h2><strong><em>Steps</em></strong></h2>  
<ol>  
    <li>Detect circRNA from RNA-seq Data.</li>  
    <li>  
        <div style="overflow-x: auto; white-space: nowrap;">  
            <code><em>sh pipeline_detection.sh -1 [sample_1.fastq or sample_1.fastq.gz] -2 [sample_2.fastq or sample_2.fastq.gz]</em></code>  
        </div>  
    </li>  
    <li>Construct a reference set of human blood full-length circRNAs.</li>  
    <li>Reconstruct incomplete circRNAs.</li>  
    <li>Annotate human blood full-length circRNAs.</li>  
    <li>Quantify human blood full-length circRNAs.</li>  
    <li>Perform Expression Analysis of human blood full-length circRNAs.</li>  
</ol>