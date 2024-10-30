<!DOCTYPE html>
<html lang="zh">

<head>
    <meta charset="UTF-8">
    <title>AQUARIUM-HB</title>
</head>

<body>

    <h1><strong><em>AQUARIUM-HB</em></strong></h1>
    <p>A bioinformatics pipeline to identify, annotate, quantify, and analyze human blood circular RNAs from RNA-seq data.</p>
    <p>For detailed documentation, please visit <a href="https://njucmbioinfo.github.io/AQUARIUM-HB/">Project Documentation</a>.</p>

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

    <h2><strong><em>Data Requirements</em></strong></h2>
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
                <code>sh AQUARIUM_HB.sh detect --fastq1 &lt;sample_1.fastq.gz&gt; --fastq2 &lt;sample_2.fastq.gz&gt; --outdir &lt;output_directory&gt;</code>
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
                <code>sh AQUARIUM_HB.sh reconstruct --cirireport_file &lt;cirireport&gt; --stoutlist_file &lt;stoutlist&gt; --reference_file &lt;reference&gt; --outputdir &lt;output_directory&gt;</code>
            </div>
        </li>
        <li>
            <strong>Quantify human blood full-length circRNAs.</strong>
            <div class="code-container">
                <code>sh AQUARIUM_HB.sh quant --fastq1 &lt;sample_1.fastq.gz&gt; --fastq2 &lt;sample_2.fastq.gz&gt; --gtfdir &lt;gtf_directory&gt; --quantdir &lt;quant_directory&gt;</code>
            </div>
        </li>
    </ul>

    <h2><strong>Usage Template</strong></h2>
    <div class="code-container">
        <code>Usage: $0 &lt;module&gt; &lt;args&gt;<br>
Modules:<br>
&nbsp;&nbsp;&nbsp;detect --fastq1 &lt;sample_1.fastq.gz&gt; --fastq2 &lt;sample_2.fastq.gz&gt; --outdir &lt;output_directory&gt;<br>
&nbsp;&nbsp;&nbsp;reference --stoutlist_path &lt;inputfile&gt; --reference_file &lt;outputfile&gt;<br>
&nbsp;&nbsp;&nbsp;reconstruct --cirireport_file &lt;cirireport&gt; --stoutlist_file &lt;stoutlist&gt; --reference_file &lt;reference&gt; --outputdir &lt;output_directory&gt;<br>
&nbsp;&nbsp;&nbsp;quant --fastq1 &lt;sample_1.fastq.gz&gt; --fastq2 &lt;sample_2.fastq.gz&gt; --gtfdir &lt;gtf_directory&gt; --quantdir &lt;quant_directory&gt;</code>
    </div>

</body>
</html>
