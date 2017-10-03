<style type="text/css">
    body{
        line-height: 2;
        font-size: 16px;
    }

    ol li{padding: 2px;}
    ul li{padding: 0px;}
    h4 {margin: 30px 0px 15px 0px;}

    div.code {
        font-family: "Courier New";
        border: 1px solid;
        border-color: #999999;
        background-color: #eeeeee;
        padding: 5px 10px;
        margin: 10px;
        border-radius: 5px;
        overflow: auto;
    }

    div.question {
        color: #666666;
        background-color: #e1eaf9;
        padding: 15px 25px;
        margin: 20px;
        font-size: 15px;
        border-radius: 20px;
    }

    div.extra {
        color: #444444;
        background-color: #e1eaf9;
        padding: 15px 25px;
        margin: 20px;
        font-size: 15px;
        border-radius: 20px;
    }

    div.question h4 {
        font-style: italic;
        margin: 10px 0px !important;
    }
</style>

<img src="../../../img/melbioinf_logo.png" height=100px>
<img src="../media/gvl_logo.jpg" height=100px align=right>

# RNA-Seq - Differential Gene Expression

**Authors: Jessica Chung, Mahtab Mirmomeni, Andrew Lonie**

-----

## Tutorial Overview

In this tutorial we cover the concepts of RNA-seq differential gene expression
(DGE) analysis using a simulated dataset from the common fruit fly, Drosophila
melanogaster.

The tutorial is designed to introduce the tools, datatypes and workflows of an
RNA-seq DGE analysis. In practice, real datasets would be much larger and
contain sequencing and alignment errors that make analysis more difficult.

In this tutorial we will:  

- introduce the types of files typically used in RNA-seq analysis
- align RNA-seq reads with an aligner (HISAT2)
- visualise RNA-seq alignment data with IGV
- use a number of different methods to find differentially expressed genes
- understand the importance of replicates for differential expression analysis

This tutorial does not cover the following steps that might do in a real
RNA-seq DGE analysis:  

- QC (quality control) of the raw sequence data
- Trimming the reads for quality and for adaptor sequences
- QC of the RNA-seq alignment data  

These steps have been omitted because the data we use in this tutorial is
synthetic and has no quality issues, unlike real data.

-----

## Learning Objectives

At the end of this tutorial you should:

 - Be familiar with basic workflow of alignment, quantification, and testing,
   for RNA-seq differential expression analysis
 - Be able to process raw RNA sequence data into a list of differentially
   expressed genes
 - Be aware of how the relationship between the number of biological replicates
   in an experiment and the statistical power available to detect differentially
   expressed genes

-----

## The data

The sequencing data you will be working with is simulated from Drosophila
melanogaster. The experiment has two conditioins, WT (wildtype) and KO
(knockout), and three samples in each condition. The sequencing data is
paired-end, so there are two files for each of the six samples. Your aim will
be to find differentially expressed genes in WT vs KO.

<img src="../media/dm_data.png" height=500px/>

-----

## Section 1: Preparation

**1.  Register as a new user in Galaxy if you don’t already have an account**

1.  Open a browser and go to a Galaxy server. This can either be your
    personal GVL server you [started previously](http://genome.edu.au/get/get#launch),
    the public [Galaxy Tutorial server](http://galaxy-tut.genome.edu.au)
    or the public [Galaxy Melbourne server](http://galaxy-mel.genome.edu.au).  
    Recommended browsers include Firefox and Chrome. Internet Explorer
    is not supported.
2.  Register as a new user by clicking **User > Register** on the top
    dark-grey bar. Alternatively, if you already have an account, login by
    clicking **User > Login**.

**2.  Import the RNA-seq data for the workshop.**

If you are using the public Galaxy Tutorial server or Galaxy Melbourne server,
you can import the data directly from Galaxy. You can do this by going to
**Shared Data > Published Histories** on the top toolbar, and selecting
the history called **RNA-Seq_Basic_2017**. Then click on "Import History" on
the top right and "start using this history" to switch to the newly imported
history.

Alternatively, if you are using your own personal Galaxy server, you can import
the data by:

1.  In the tool panel located on the left, under Basic Tools select **Get
    Data > Upload File**. Click on the **Paste/Fetch data** button on the
    bottom section of the pop-up window.
2.  Upload the sequence data by pasting the following links into the text
    input area.
    These six files are three paired-end samples from the WT flies. Make sure
    the type is specified as 'fastqsanger' when uploading.
    <div class="code">
    https://swift.rc.nectar.org.au:8888/v1/AUTH_377/public/rna_seq_basic/WT_01_R1.fastq
    <br>
    https://swift.rc.nectar.org.au:8888/v1/AUTH_377/public/rna_seq_basic/WT_01_R2.fastq
    <br>
    https://swift.rc.nectar.org.au:8888/v1/AUTH_377/public/rna_seq_basic/WT_02_R1.fastq
    <br>
    https://swift.rc.nectar.org.au:8888/v1/AUTH_377/public/rna_seq_basic/WT_02_R2.fastq
    <br>
    https://swift.rc.nectar.org.au:8888/v1/AUTH_377/public/rna_seq_basic/WT_03_R1.fastq
    <br>
    https://swift.rc.nectar.org.au:8888/v1/AUTH_377/public/rna_seq_basic/WT_03_R2.fastq
    <br>
    </div>
    These six files are three paired-end samples from the KO flies.
    Make sure the type is specified as 'fastqsanger' when uploading.
    <div class="code">
    https://swift.rc.nectar.org.au:8888/v1/AUTH_377/public/rna_seq_basic/KO_01_R1.fastq
    <br>
    https://swift.rc.nectar.org.au:8888/v1/AUTH_377/public/rna_seq_basic/KO_01_R2.fastq
    <br>
    https://swift.rc.nectar.org.au:8888/v1/AUTH_377/public/rna_seq_basic/KO_02_R1.fastq
    <br>
    https://swift.rc.nectar.org.au:8888/v1/AUTH_377/public/rna_seq_basic/KO_02_R2.fastq
    <br>
    https://swift.rc.nectar.org.au:8888/v1/AUTH_377/public/rna_seq_basic/KO_03_R1.fastq
    <br>
    https://swift.rc.nectar.org.au:8888/v1/AUTH_377/public/rna_seq_basic/KO_03_R2.fastq
    <br>
    </div>
    Then, upload this file of gene definitions. You don't need to specify
    the type for this file as Galaxy will auto-detect the file as a GTF
    file.
    <div class="code">
    https://swift.rc.nectar.org.au:8888/v1/AUTH_377/public/rna_seq_basic/ensembl_dm3.chr4.gtf
    </div>
    You should now have 13 files in your history.

    **Note:** If you log out of Galaxy and log back at a later time your data
    and results from previous experiments will be available in the right panel
    of your screen called the ‘History’


**3.  View and have an understanding of the files involved in RNA-seq analysis.**

1.  You should now have the following files in your Galaxy history:

    6 files containing paired-ended reads for the WT samples:
    <ul>
    <li>WT_01_R1.fastq  
    <li>WT_01_R2.fastq  
    <li>WT_02_R1.fastq  
    <li>WT_02_R2.fastq  
    <li>WT_03_R1.fastq  
    <li>WT_03_R2.fastq  
    </ul>

    6 files containing paired-ended reads for the KO samples:
    <ul>
    <li>KO_01_R1.fastq  
    <li>KO_01_R2.fastq  
    <li>KO_02_R1.fastq  
    <li>KO_02_R2.fastq  
    <li>KO_03_R1.fastq  
    <li>KO_03_R2.fastq  
    </ul>

    And 1 gene annotation file:
    <ul>
    <li>ensembl_dm3.chr4.gtf  
    </ul>

    These files can be renamed by clicking the **pen icon** if you wish.

2.  These 12 sequencing files are in FASTQ format and have the file
    extension: .fastq. If you are not familiar with the FASTQ format, [click
    here for an overview](https://en.wikipedia.org/wiki/FASTQ_format).  

    Each condition has three samples, and each sample has two files (an R1
    file containing forward reads and an R2 file containing reverse reads).

    Click on the **eye icon** to the top right of any FASTQ file to view the
    first part of the file.

    **Note:** Since the reads in this dataset are synthetic, they do not have
    real quality scores.

    **Note:** The reads are paired-end, i.e. WT_01_R1.fastq and
    WT_01_R2.fastq are paired reads from one sequencing run. If you're
    unfamiliar with paired-end sequencing, you can read about it
    [here](https://www.illumina.com/science/technology/next-generation-sequencing/paired-end-vs-single-read-sequencing.html).

3.  The gene annotation file (ensembl_dm3.chr4.gtf) is in GTF format. This file
    describes where the genes are located in the D. melanogaster reference
    genome, filtered for genes on chromosome 4.
    Each feature is defined by a chromosomal start and end point, feature type
    (CDS, gene, exon etc), and parent gene and transcript.
    We will examine this file more closely later in Section 3 of this
    tutorial.
    More information on the GTF format can be found
    [here](http://asia.ensembl.org/info/website/upload/gff.html).

-----

## Section 2: Alignment with HISAT2

In this section we map the reads in our FASTQ files to a reference genome. As
these reads originate from mRNA, we expect some of them will cross exon/intron
boundaries when we align them to the reference genome. We will use HISAT to
perform our alignment. HISAT2 is a fast, splice-aware, alignment program that
is a successor to TopHat2. More information on
HISAT2 can be found [here](https://ccb.jhu.edu/software/hisat2/index.shtml).

**1.  Align the RNA-seq reads to a reference genome.**

In the left tool panel menu, under NGS Analysis, select
**NGS: RNA Analysis > HISAT2** and set the parameters as follows:  

- **Input data format** FASTQ
- **Single end or paired reads?** Individual paired reads
- **Forward reads:**  
(Click on the **multiple datasets icon** and select all six of the forward
FASTQ files ending in \*1.fastq. This should be correspond to every
second file (1,3,5,7,9,11). This can be done by holding down the
ctrl key (Windows) or the command key (OSX) to select multiple files.)
    - WT_01_R1.fastq
    - WT_02_R1.fastq
    - WT_03_R1.fastq
    - KO_01_R1.fastq
    - KO_02_R1.fastq
    - KO_03_R1.fastq

- **Reverse reads:**  
(Click on the **multiple datasets icon** and select all six of the reverse
FASTQ files ending in \*2.fastq.)  
    - WT_01_R2.fastq
    - WT_02_R2.fastq
    - WT_03_R2.fastq
    - KO_01_R2.fastq
    - KO_02_R2.fastq
    - KO_03_R2.fastq
- **Source for the reference genome to align against:** Use
built-in genome
- **Select a reference genome:** D. melanogaster Apr. 2006 (BDGP R5/dm3)
  (dm3)
- Use defaults for the other fields
- Execute

<img src="../media/rna_basic_hisat2.png"/>

Note: This may take a few minutes, depending on how busy the server is.

**2.  Examine the alignment stats**

HISAT2 outputs one bam file for each set of paired-end read files. Rename the 6
files into a more meaningful name (e.g. 'HISAT on data 2 and data 1' to 'WT_01.bam')
by using the **pen icon** next to the file.

These files are BAM files (short for
[Binary Alignment Map](https://en.wikipedia.org/wiki/Binary_Alignment_Map))
and like the name suggests, is a binary file. This means we can't use the
eye icon to view the data in Galaxy; we need to use software that can read the
file or convert it into it's plain-text equivalent (SAM) to view it as text.
In section 3, we'll use a genome viewer to view our alignments.

HISAT2 also outputs some information to stderr which we can preview by
clicking on the dataset name. To view the raw file, click the "info" button
(view details) of a dataset, say WT_01.bam, and find the "Tool Standard Error"
row under "Job Information" in the table. Click the "stderr" link to view
the alignment summary output.

```
16046 reads; of these:
  16046 (100.00%) were paired; of these:
    104 (0.65%) aligned concordantly 0 times
    13558 (84.49%) aligned concordantly exactly 1 time
    2384 (14.86%) aligned concordantly >1 times
    ----
    104 pairs aligned concordantly 0 times; of these:
      1 (0.96%) aligned discordantly 1 time
    ----
    103 pairs aligned 0 times concordantly or discordantly; of these:
      206 mates make up the pairs; of these:
        106 (51.46%) aligned 0 times
        91 (44.17%) aligned exactly 1 time
        9 (4.37%) aligned >1 times
99.67% overall alignment rate
```

Here we see we have a very high alignment rate, which is expected since the
data we have is simulated and has no contamination.

-----

## Section 3: Visualise the aligned reads

The purpose of this step is to :

- visualise the quantitative, exon-based nature of RNA-seq data
- visualise the expression differences between samples represented by the
  quantity of reads, and
- become familiar with the [Integrative Genomics Viewer
  (IGV)](https://www.broadinstitute.org/igv/)-- an interactive
  visualisation tool by the Broad Institute.  


To visualise the alignment data:

1.  Click on one of the BAM files, for example 'WT_01.bam'.
2.  Click on Display with IGV **'webcurrent'** (or 'local' if you have IGV
    installed on your computer. You will need to open IGV before you click on
    'local'). This should download a .jnlp Java Web Start file to your
    computer. Open this file to run IGV. (You will need Java installed on your
    computer to run IGV)
3.  Once IGV opens, it will show you the BAM file. (Note:
    this may take a bit of time as the data is downloaded to IGV)
4.  Select **chr4** from the second drop box under the toolbar. Zoom in to
    view alignments of reads to the reference genome.
    You should see the characteristic distribution of RNA-seq reads across
    the exons of the genes, with some gaps at intron/exon boundaries.
    The number of reads aligned to a particular gene is proportional to the
    abundance of the RNA derived from that gene in the sequenced sample.
    (Note that IGV already has a list of known genes of most major organisms
    including Drosophila, which is why you can see the genes in the bottom
    panel of IGV.)
5.  View differentially expressed genes by viewing two alignment files
    simultaneously. The aim of this tutorial is to statistically test
    differential expression, but first it’s useful to reassure ourselves
    that the data looks right at this stage by comparing the aligned reads
    for condition 1 (WT) and condition 2 (KO).

    Select 'KO_02.bam' and click on 'display with IGV local'. This time we are
    using the **'local'** link, as we already have an IGV window up and running
    locally from the last step. Once the file has loaded, try to find some
    genes that look differentially expressed.

    <img src="../media/rna_igv.png" height=400px/>

    If you can't find any, try changing the location to
    **chr4:816349-830862** using the field on the top toolbar.
    The 'Sox102F' gene in this area looks like it has many more reads
    mapped in WT than in KO. Hover over the coverage track to view the read
    depth of the area. But, of course, it may be that there are many more
    reads in the library for WT than KO. So we need to statistically
    normalise the read counts before we can say anything definitive,
    which we will do in the next section.

**[Optional]** Visualise the aligned reads in Trackster  
If you have trouble getting IGV to work, you can also use the inbuilt
Galaxy genome browser, Trackster, to visualise alignments. Trackster has fewer
features than IGV, but sometimes it may be more convenient to use as it only
requires the browser.

1.  On the top bar of Galaxy, select **Visualization > New Track Browser**.
2.  Name your new visualization and select D. melanogaster (dm3) as the
    reference genome build.
3.  Click the **Add Datasets to Visualization** button and select WT_01.bam and
    KO_01.bam by using the checkboxes on the left.
4.  Select chr4 from the dropdown box. You can zoom in and out using the
    buttons on the top toolbar. You can also add more tracks using the Add
    Tracks icon located on the top right.
5.  Next to the drop down list, click on the chromosomal position number
    display and specify the location **chr4:816349-830862**.  

Before starting the next section, leave the Trackster interface and return
to the analysis view of Galaxy by clicking 'Analyze Data' on the top
Galaxy toolbar.

-----


## Section 4. Quantification

HTSeq-count counts the number of the reads from each bam file that map to the
genomic features in the provided annotation file. For each feature (a
gene for example) we will obtain a numerical value associated with the
expression of that feature in our sample (i.e. the number of reads that
were aligned to that gene).

**1.  Examine the GTF file**

Click on the **eye icon** to display the ensembl_dm3.chr4.gtf file in Galaxy.

This GTF file is essentially a list of chromosomal features
which together define genes. Each feature is in turn defined by a
chromosomal start and end point, feature type (CDS, gene, exon etc),
and parent gene and transcript. Importantly, a gene may have many features,
but one feature will belong to only one gene.
More information on the GTF format can be found
[here](http://asia.ensembl.org/info/website/upload/gff.html).
<br><img src="../media/rna_basic_gtf.png" height=200px style="margin: 10px;"/><br>
The ensembl_dm3.chr4.gtf file contains ~4900 features which together define
the 92 known genes on chromosome 4 of Drosophila melanogaster.

**2.  Run HTSeq-count**

1.  Use HTSeq-count to count the number of reads for each feature.  
    In the left tool panel menu, under NGS Analysis, select
    **NGS: RNA Analysis > htseq-count** and set the parameters as follows:  
    - **Aligned SAM/BAM File:**  
      (Select 'Multiple datasets', then select all six bam files using the shift key.)
        - WT_01.bam
        - WT_02.bam
        - WT_03.bam
        - KO_01.bam
        - KO_02.bam
        - KO_03.bam
    - **GFF File:** ensembl_dm3.chr4.gtf
    - **Stranded:** No
    - **ID Attribute:** gene_name
    - Use defaults for the other fields
    - Execute

2.  In the previous step, each input BAM file outputted two files. The first
    file contains the counts for each of our genes. The second file
    (ending with "(no feature)") contains the stats for the reads that weren't
    able to be uniquely aligned to a gene. We don't need the "(no feature)"
    files so we can remove then with the delete "X" button on the top right.

3.  Rename the remaining six files from htseq-count to meaningful names,
    such as WT_01, WT_02, etc.

**3.  Generate a count matrix**

1.  Generate a combined count matrix by combining our six files.
    In the left tool panel menu, under NGS Analysis, select
    **NGS: RNA Analysis > Generate count matrix** and set the parameters as follows:  
    - **Count files from your history:**  
        (Select all six count files using the shift key.)
        - WT_01
        - WT_02
        - WT_03
        - KO_01
        - KO_02
        - KO_03
    - Use defaults for the other fields
    - Execute

Examine the outputted matrix by using the **eye icon**.  
Each column corresponds to a sample and each row corresponds to a gene. By
sight, see if you can find a gene you think is differentially expressed
just by looking at the counts.

We now have a count matrix which we will now use to find differentially
expressed genes between WT samples and KO samples.

-----

## Section 5. Degust

Degust is an interactive visualiser for analysing RNA-seq data. It runs as a
web service and can be found at [degust.erc.monash.edu/](http://degust.erc.monash.edu/).

<img src="../media/degust_01.png" height=400px style="display:block; margin-left: auto; margin-right:auto;">

**1. Load count data into Degust**

1.  In Galaxy, download the count matrix you generated in the last section
    using the **disk icon**.
2.  Go to [degust.erc.monash.edu/](http://degust.erc.monash.edu/)
    and click on "Upload your counts file".
3.  Click "Choose file" and upload the recently downloaded Galaxy tabular file
    containing your RNA-seq counts.

**2. Configure your uploaded data**

1.  Give your visualisation a name.
2.  For the Info column, select "gene_id".
3.  Add two conditions: WT and KO. For each condition, select the three
    samples which correspond with the condition.
4.  Set min gene CPM to 1 in at least 3 samples.
4.  Click **Save changes** and view your data.

<img src="../media/degust_02.png" style="display:block; margin-left: auto; margin-right:auto;">

\endshowable


Read through the Degust tour of features. Explore the parallel coordinates plot,
MA plot, MDS plot, heatmap and gene list. Each is fully interactive and
influences other portions on the display depending on what is selected.

<img src="../media/degust_03.png" style="display:block; margin-left: auto; margin-right:auto;">

On the right side of the page is an options module which can set thresholds to
filter genes using statistical significance or absolute-fold-change.

On the left side is a dropdown box you can specify the method (Voom/Limma or
edgeR) used to perform differential expression analysis on the data. You can
also view the R code by clicking "Show R code" under the options module on
the right.

**4. Explore the demo data**

Degust also provides an example dataset with 4 conditions and more genes. You
can play with the demo dataset by clicking on the "Try the demo" button on the
Degust homepage. The demo dataset includes a column with an EC number for each
gene. This means genes can be displayed on Kegg pathways using the module on
the right.


-----

## Section 5. DESeq2

In this section we'll use the "DESeq2" tool in Galaxy to do our differential
gene analysis. This tool uses the separate HTSeq files we generated in section
4.

Similar to Voom/Limma or edgeR that was used in Degust to statistically test
our data, DESeq2 will:

  - statistically test for expression differences in normalised read counts for
    each gene, taking into account the variance observed between samples,  
  - for each gene, calculate the p-value of the gene being differentially
    expressed-- this is the probability of seeing the data or something more
    extreme given the null hypothesis (that the gene is not differentially
    expressed between the two conditions),
  - for each gene, estimate the fold change in expression between the two
    conditions.

<p></p>

1.  Use DESeq2 to find differentially expressed features from the count data.  
    In the left tool panel menu, under NGS Analysis, select **NGS: RNA Analysis
    > DESeq2** and set the parameters as follows:
    - **1: Factor**
        - **Specify a factor name:** condition
        - **1: Factor level:**
            - **Specify a factor level:** WT  
            (Select the three WT htseq-count files.)
                - WT_01
                - WT_02
                - WT_03
        - **2: Factor level:**
            - **Specify a factor level:** KO  
            (Select the three KO htseq-count files.)
                - KO_01
                - KO_02
                - KO_03
        - Use defaults for the other fields
        - Execute

<img src="../media/rna_basic_deseq2.png"/>

2.  Have a look at the outputs of DESeq2. We will now filter significant
    (adjusted p-value < 0.05) genes from the DESeq2 result file.  
    Under Basic Tools, click on **Filter and Sort > Filter**:
    - **Filter:** "DESeq2 result file on data ..."
    - **With following condition:** c7 < 0.05
    - Execute

How many differentially expressed genes with adjusted p-value < 0.05 are there?

-----

## Section 6. The importance of replicates

1.  Repeat the previous differential expression analysis with two samples in
    each group instead of three. How do you expect your results to differ when
    using fewer samples?

2.  Filter genes with adjusted-p-value < 0.05. How many genes are significant?

3.  Run DESeq2 again, using only one sample from each group. How many genes are
    now significant?

4.  Can you find genes that were identified as differentially expressed when
    using three samples in each condition that were not identified as
    differentially expressed when using two samples? What do you expect these
    gene's counts or logFC values to look like compared to genes that remained
    statistically significance? Have a look at the counts or the logFC values
    of these genes.


The identification of differentially expressed genes is based on the size
of the difference in expression and the variance observed across multiple
replicates. This demonstrates how important it is to have biological replicates
in differential gene expression experiments.

Myoglianin is an example of a gene that showed up as differentially expressed
when we did a 3 vs 3 comparsion but not with a 2 vs 2 comparsion.
If we say that genes like Myoglianin was **truly** differentially expressed, we
can call these instances where the true differentially expressed genes are not
identified as false negatives. Generally, increasing replicates decreases the
number of false negatives.

It is also more likely to see more false positives when using an insufficient
number of replicates. False positives can be defined as identifiying a gene as
differentially expressed when it is, in reality, not.

-----

## Optional extension

Have a go at doing another differential expression analysis with the following
Saccharomyces cerevisiae data from chromosome I. This time, the two conditions
are called 'batch' and 'chem', and like before, there are three samples per
condition.

Batch sequence data:
<div class="code">
https://swift.rc.nectar.org.au:8888/v1/AUTH_a3929895f9e94089ad042c9900e1ee82/RNAseqDGE_ADVNCD/batch1_chrI_1.fastq
<br>
https://swift.rc.nectar.org.au:8888/v1/AUTH_a3929895f9e94089ad042c9900e1ee82/RNAseqDGE_ADVNCD/batch1_chrI_2.fastq
<br>
https://swift.rc.nectar.org.au:8888/v1/AUTH_a3929895f9e94089ad042c9900e1ee82/RNAseqDGE_ADVNCD/batch2_chrI_1.fastq
<br>
https://swift.rc.nectar.org.au:8888/v1/AUTH_a3929895f9e94089ad042c9900e1ee82/RNAseqDGE_ADVNCD/batch2_chrI_2.fastq
<br>
https://swift.rc.nectar.org.au:8888/v1/AUTH_a3929895f9e94089ad042c9900e1ee82/RNAseqDGE_ADVNCD/batch3_chrI_1.fastq
<br>
https://swift.rc.nectar.org.au:8888/v1/AUTH_a3929895f9e94089ad042c9900e1ee82/RNAseqDGE_ADVNCD/batch3_chrI_2.fastq
<br>
</div>
Chem sequence data:
<div class="code">
https://swift.rc.nectar.org.au:8888/v1/AUTH_a3929895f9e94089ad042c9900e1ee82/RNAseqDGE_ADVNCD/chem1_chrI_1.fastq
<br>
https://swift.rc.nectar.org.au:8888/v1/AUTH_a3929895f9e94089ad042c9900e1ee82/RNAseqDGE_ADVNCD/chem1_chrI_2.fastq
<br>
https://swift.rc.nectar.org.au:8888/v1/AUTH_a3929895f9e94089ad042c9900e1ee82/RNAseqDGE_ADVNCD/chem2_chrI_1.fastq
<br>
https://swift.rc.nectar.org.au:8888/v1/AUTH_a3929895f9e94089ad042c9900e1ee82/RNAseqDGE_ADVNCD/chem2_chrI_2.fastq
<br>
https://swift.rc.nectar.org.au:8888/v1/AUTH_a3929895f9e94089ad042c9900e1ee82/RNAseqDGE_ADVNCD/chem3_chrI_1.fastq
<br>
https://swift.rc.nectar.org.au:8888/v1/AUTH_a3929895f9e94089ad042c9900e1ee82/RNAseqDGE_ADVNCD/chem3_chrI_2.fastq
<br>
</div>
Gene annotation:
<div class="code">
https://swift.rc.nectar.org.au:8888/v1/AUTH_a3929895f9e94089ad042c9900e1ee82/RNAseqDGE_ADVNCD/genes.gtf
</div>
