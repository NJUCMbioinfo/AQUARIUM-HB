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

<h3><strong>Required R Packages</strong></h3>  
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

<h2><strong>Steps</strong></h2>  
<ul>  
    <li>  
        <strong>Detect circRNA from RNA-seq Data.</strong>  
        <div class="code-container">  
            <code>sh AQUARIUM_HB.sh detect --fastq1 &lt;sample_1.fastq.gz&gt; --fastq2 &lt;sample_2.fastq.gz&gt; --outdir &lt;outdir&gt;</code>  
        </div>  
    </li>  
    <li>  
        <strong>Construct a reference set of human blood full-length circRNAs.</strong>  
        <div class="code-container">  
            <code>sh AQUARIUM_HB.sh reference --stoutlist_path &lt;inputfile&gt; --reference_file &lt;outputfile&gt;</code>  
        </div>  
    </li>  
    <li>  
        <strong>Reconstruct incomplete circRNAs.</strong>  
        <div class="code-container">  
            <code>sh AQUARIUM_HB.sh reconstruct --cirireport_file &lt;cirireport&gt; --stoutlist_file &lt;stoutlist&gt; --reference_file &lt;reference&gt; --outputdir &lt;outputdir&gt;</code>  
        </div>  
    </li>  
    <li>  
        <strong>Quantify human blood full-length circRNAs.</strong>  
        <div class="code-container">  
            <code>sh AQUARIUM_HB.sh quant --fastq1 &lt;sample_1.fastq.gz&gt; --fastq2 &lt;sample_2.fastq.gz&gt; --gtfdir &lt;gtfdir&gt; --quantdir &lt;quantdir&gt;</code>  
        </div>  
    </li>  
</ul>  

</body>  
</html>