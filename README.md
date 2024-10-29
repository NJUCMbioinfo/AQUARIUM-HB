<!DOCTYPE html>  
<html lang="zh">  
<head>  
    <meta charset="UTF-8">  
    <meta name="viewport" content="width=device-width, initial-scale=1.0">  
    <title>AQUARIUM-HB</title>  
    <style>  
        ul, ol {  
            list-style-type: none; /* 去掉默认的项目符号或编号 */  
            padding-left: 0; /* 去掉左侧缩进 */  
        }  
    </style>  
</head>  
<body>  

<h1><strong><em>AQUARIUM-HB</em></strong></h1>  

<p>A bioinformatics pipeline to identify, annotate, quantify, and analyze human blood circular RNAs from RNA-seq data.</p>  

<p>For detailed documentation, please visit the <a href="https://njucmbioinfo.github.io/AQUARIUM-HB/">Project Documentation</a>.</p>  

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
<ul>  
    <li>Detect circRNA from RNA-seq Data.</li>  
    <li>  
        <div style="overflow-x: auto; white-space: nowrap; width: 100%; border: 1px solid #ccc; padding: 5px;">  
            <code><em>sh pipeline_detection.sh -1 sample_1.fastq -2 sample_2.fastq</em></code>    
        </div>  
    </li>  
    <li>Construct a reference set of human blood full-length circRNAs.</li>  
    <li>Reconstruct incomplete circRNAs.</li>  
    <li>Annotate human blood full-length circRNAs.</li>  
    <li>Quantify human blood full-length circRNAs.</li>  
    <li>Perform Expression Analysis of human blood full-length circRNAs.</li>  
</ul>  

</body>  
</html>