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

<h2><strong><em>Steps</em></strong></h2>  
<ul>  
    <li>Detect circRNA from RNA-seq Data.</li>  
    <li>  
        <div style="overflow-x: auto; white-space: nowrap; width: 100%; border: 1px solid #ccc; padding: 5px;">  
            <code><em>sh AQUARIUM_HB.sh -1 sample_1.fastq -2 sample_2.fastq</em></code>  
        </div>  
    </li>  
    <li>Construct a reference set of human blood full-length circRNAs.</li>  
    <li>  
        <div style="overflow-x: auto; white-space: nowrap; width: 100%; border: 1px solid #ccc; padding: 5px;">  
            <code><em>sh AQUARIUM_HB.sh reference -i stout.list_path.txt -o ReferenceSet.txt</em></code>  
        </div>  
    </li>  
    <li>Reconstruct incomplete circRNAs.</li>  
    <li>  
        <div style="overflow-x: auto; white-space: nowrap; width: 100%; border: 1px solid #ccc; padding: 5px;">  
            <code><em>sh AQUARIUM_HB.sh reconstruct -r ReferenceSet.txt -s sampleA</em></code>  
        </div>  
    </li>  
    <li>Quantify human blood full-length circRNAs.</li>  
    <li>  
        <div style="overflow-x: auto; white-space: nowrap; width: 100%; border: 1px solid #ccc; padding: 5px;">  
            <code><em>sh AQUARIUM_HB.sh quant -s sampleA</em></code>  
        </div>  
    </li>  
</ul>  

</body>  
</html>
