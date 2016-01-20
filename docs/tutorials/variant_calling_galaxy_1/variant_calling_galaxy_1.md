<p>
<img src="../media/vlsci_logo.jpg" alt="VLSCI logo" align="left" width="164"/>
<img src="../media/gvl_logo.jpg" alt="GVL logo" align="right" width="112"/>
</p>
<p><br></p>
<p><br></p>

# Table of Contents
1. [Tutorial Overview](#1-tutorial-overview)
2. [Background](#2-background)
3. [Preparation](#3-preparation)
4. [Quality Control](#4-quality-control)
5. [Alignment to the reference - (FASTQ to BAM)](#5-alignment-to-the-reference---fastq-to-bam))
6. [Calling single nucleotide variations(SNVs)](#6-calling-single-nucleotide-variations-snvs)
7. [Calling small insertions and deletetions](#7-calling-small-insertions-and-deletions-indels)


## 1. Tutorial Overview

In this tutorial we cover the concepts of detecting small variants (SNVs and indels) in human genomic DNA using a small set of reads from chromosome 22.

Note:

The tutorial is designed to introduce the tools, datatypes and workflow of variation detection. We filter the variations manually to understand what is actually happening in variant calling. In practice the datasets would be much larger and you would use more sophisticated tools to call, annotate and filter variants.

## 2. Background

Some background reading material - [background]

#### Where is the data in this tutorial from?

The workshop is based on analysis of short read data from the exome of chromosome 22 of a single human individual. There are one million 76bp reads in the dataset, produced on an Illumina GAIIx from exome-enriched DNA. This data was generated as part of the [1000 genomes] Genomes project.

## 3. Preparation

1. **Make sure you have an instance of Galaxy ready to go.**
    *  If not - go to our [Melbourne Galaxy] instance.
2. **Import data for the tutorial.**
    * In this case, we are uploading a [FASTQ] file.
    * **Method 1**  
        *Paste/Fetch data from a URL to Galaxy.*
        1. In the Galaxy tools panel (left), click on *Get Data* and choose *Upload File*.
        2. Click *Paste/Fetch data* and paste the following URL into the box
<div class="code">
https://swift.rc.nectar.org.au:8888/v1/AUTH_a3929895f9e94089ad042c9900e1ee82/VariantDet_BASIC/NA12878.GAIIx.exome_chr22.1E6reads.76bp.fastq
<br>
</div>
        3. Select *Type* as **fastqsanger** and click *Start*.
        4. Once the upload status turns *green*, it means the upload is complete. You should now be able to see the file in the Galaxy history panel (right).

    Alternatively, if you have a local file to upload (*For the purpose of this tutorial we can stick with the option above*):

    * **Method 2**  
        *Upload data to Galaxy.*
        1. In the Galaxy tools panel (left), click on *Get Data* and choose *Upload File*.
        2. From *Choose local file* select the downloaded FASTQ file. Select *Type* as **fastqsanger** and click *Start*.
        3. Once the upload status turns *green*, it means the upload is complete. You should now be able to see the file in the Galaxy history panel (right).


Summary:  
So far, we have started a Galaxy instance, got hold of our data and uploaded it to
the Galaxy instance.  
Now we are ready to perform our analysis.


## 4. Quality Control

The first step is to evaluate the quality of the raw sequence data. If the quality is poor, then adjustments can be made - e.g. trimming the short reads, or adjusting your expectations of the final outcome!

#### Take a look at a FASTQ file
1. Click on the eye icon to the top right of the fastq file to view the a snippet of the file.
2. Note that each read is represented by 4 lines:
  * read identifier
  * short read sequence
  * separator
  * short read sequence quality scores

```text
e.g.
identifier:    @61CC3AAXX100125:7:72:14903:20386/1
read sequence: TTCCTCCTGAGGCCCCACCCACTATACATCATCCCTTCATGGTGAGGGAGACTTCAGCCCTCAATGCCACCTTCAT
separator:     +
quality score: ?ACDDEFFHBCHHHHHFHGGCHHDFDIFFIFFIIIIHIGFIIFIEEIIEFEIIHIGFIIIIIGHCIIIFIID?@<6

```

For more details see [FASTQ].

#### Assessing read quality from the FASTQ files
From the Galaxy tools panel, select
```text

NGS: QC and manipulation > FastQC: Comprehensive QC

The input FASTQ file will be selected by default. Keep the other defaults and click 'Execute'

```

Note the batch processing interface of Galaxy:

 * grey = waiting in queue
 * yellow = running
 * green = finished
 * red = tried to run and failed

Wait for job to be complete.

Click on the eye icon to view the newly generated data (in this case a set of quality metrics for the FASTQ data).
Look at the various quality scores. The data looks pretty good - *high Per base sequence quality* (avg. above 30).

*Note that the 'Sequence Duplication Levels' are marked as high. Normally we would run another tool to remove duplicates (technical PCR artifacts) but for the sake of brevity we will omit this step.*

## 5. Alignment to the reference - (FASTQ to BAM)

The basic process here to map individual reads - from the input sample FASTQ file - to a matching region on the reference genome.

1.  Map/align the reads with the [BWA] tool to Human reference genome 19 (hg19) [UCSC hg19].
    From the Galaxy tools panel, select
<div class=code>
NGS: Mapping > Map with BWA-MEM [3-5mins]<br>
<br>
From the options:<br>
Using reference genome: set to hg19.<br>
Single or Paired-end reads: set to Single<br>
<br>
keep other options as default and click execute<br>
</div>

    *Note: This is the longest step in the workshop and will take a few minutes, possibly more depending on how many people are also scheduling mappings*

2.  Sort the BAM file.
  From the Galaxy tools panel, select
<div class=code>
NGS: SAM Tools > Sort BAM dataset
<br><br>
From the options:
<br>
BAM File: set to the output from the alignment BAM file
<br>
Sort by: Chromosomal coordinates
<br><br>
Keep other options as default and click execute
<br>
</div>

3.  To examine the output sorted BAM file, we need to first convert it into readable [SAM] format.
  From the Galaxy tools panel, select
<div class=code>
NGS: SAM Tools > BAM-to-SAM
<br><br>
From the options:
<br>
BAM File to Convert: set to the output to the sorted BAM file
<br><br>
Keep other options as default and click execute
<br>
</div>

4.  Examine the generated Sequence Alignment Map (SAM) file.
  * Click the eye icon next to the newly generated file
  * Familiarise yourself with the [SAM] format
  * Note that some reads have mapped to non-chr22 chromosomes (see column 3).

    This is the essence of alignment algorithms - the aligner does the best it can, but because of compromises in accuracy vs performance and repetitive sequences in the genome, not all the reads will necessarily align to the ‘correct’ sequence or could this be suggesting the presence of a structural variant?

    >  <img src="../media/tips.png" alt="Tip" height="42" width="42"/>
    Tip:
    Galaxy auto-generates a name for all outputs. Therefore, it is advisable to choose a more meaningful name to these outputs.
    >
    > This can be done as follows:
    >
    > Click on the pencil icon (edit attributes) and change Name e.g. Sample.bam or Sample.sam or Sample.sorted.bam etc.


5.  Assessing the alignment data

    *Mapping statistics*

    From the Galaxy tools panel, select
<div class=code>
NGS: SAM Tools > IdxStats
<br><br>
From the options:
<br>
The BAM: select the sorted BAM file
<br><br>
keep other options as default and click execute
<br>
</div>
    The output format (tab delimited):

    1. reference sequence identifier (chromosome)
    2. reference sequence length
    3. number of mapped reads
    4. number of placed but unmapped reads  

    *Check the number of aligned and unaligned reads + other quality metrics*

    From the Galaxy tools panel, select
<div class=code>
NGS: Sam Tools > Flagstat
<br><br>
From the options:
<br>
The BAM: select the sorted BAM file
<br><br>
keep other options as default and click execute
<br>
</div>

    Note that in this case the statistics are not very informative. This is because the dataset has been generated for this workshop and much of the noise has been removed (and in fact we just removed a lot more noise in the previous step); also we are using single ended read data rather than paired-end so some of the metrics are not relevant.

6.  Visualize the BAM file.
    * Download the sorted BAM file
    From the Galaxy history panel, select
<div class=code>
Click on the sorted BAM file.<br>
Click on the disk icon > Download dataset<br>
Click on the disk icon > Download bam_index<br>
<br>
Result:<br>
This will result in two files being downloaded 1) bam file 2) bam index file<br>
<br>
Open [IGV] browser > Open downloaded BAM file > select chromosome 22<br>
<br>
Can't see anything!<br>
<br>
Zoom in further or type a gene of interest and now you should see aligned reads.<br>
</div>

<a href="../media/igv1.jpg"><img src="../media/igv1.jpg" alt="IGV view 1" width="640px" style="display:block; margin-left: auto; margin-right:auto;"/></a>

Try looking at region chr22:36,006,744-36,007,406
Can you see a few variants?  

<a href="../media/igv_mb.jpg"><img src="../media/igv_mb.jpg" alt="IGV view 1" width="640px" style="display:block; margin-left: auto; margin-right:auto;"/></a>


## 6. Calling single nucleotide variations (SNVs)

1.  Generate a pileup file from the aligned reads (sorted bam file previous step 2). A pileup is essentially a column wise representation of the aligned read - at the base level - to the reference. The pileup file summarises all data from the reads at each genomic region that is covered by at least one read.

    From the Galaxy tools panel, select
<div class=code>
NGS: SAMtools>Generate Pileup<br>
<br>
From the options:<br>
Call consensus according to MAQ model = Yes<br>
This generates a called 'consensus base' for each chromosomal position.<br>
<br>
keep other options as default and click execute<br>
</div>

    For each output file, Galaxy tries to assign a datatype attribute to every file. For each output file, Galaxy tries to assign a datatype attribute to every file.

    Firstly, remember to rename output files to something more meaningful e.g. the workflow stage

    1. For the above output, click on the pencil icon (edit attributes)
    2. Choose Datatype link from the top column.
    3. Although it is a tabular file, for downstream processing we want to tell Galaxy that this is a *pileup* file.
    4. From the drop-down, select Pileup and click save.

    >  <img src="../media/tips.png" alt="Tip" height="42" width="42"/>
    >  Tip:
    >Get familiar with the [Pileup] format.<br>
    >* 1. chromosome<br>
    >* 2. position<br>
    >* 3. current reference base<br>
    >* 4. consensus base from the mapped reads<br>
    >* 5. consensus quality<br>
    >* 6. SNV quality<br>
    >* 7. maximum mapping quality<br>
    >* 8. coverage<br>
    >* 9. bases within reads<br>
    >* 10. quality values<br>
    >Further on (9):<br>
    >   * Each character represents one of the following: the longer this string, higher the coverage<br>
    >   * . = match on forward strand for that base<br>
    >   * , = match on reverse strand<br>
    >   * ACGTN = mismatch on forward<br>
    >   * acgtn = mismatch on reverse<br>
    >   * +[0-9]+[ACGTNacgtn]+' = insertion between this reference position and the next<br>
    >   * -[0-9]+[ACGTNacgtn]+' = deletion between this reference position and the next<br>
    >   * ^ = start of read<br>
    >   * $ = end of read<br>
    >   * BaseQualities = one character per base in ReadBases, ASCII encoded Phred scores<br>

2.  Filtering pileup to get a list of SNVs

    The pileup from the above steps described *all* aligned positions. For the purpose of calling SNVs, we are after:

    * aligned positions where the consensus base is difference from the reference base
    * the different consensus base is covered by at least e.g. 10 reads    
    * removing poor quality bases

    From the Galaxy tools panel, select
<div class=code>
NGS: SAM Tools > Filter Pileup<br>
<br>
From the options:<br>
which contains = Pileup with ten columns (with consensus)<br>
Do not report positions with coverage lower than = 10<br>
Convert coordinates to intervals = Yes<br>
<br>
keep other options as default and click execute<br>
</div>

    How many variants (*intervals*) do you see?

    If you see ~30,000 variants - do not be surprised. Our filtering criteria was not very stringent.

    > <img src="../media/tips.png" alt="Tip" height="42" width="42"/>
    > Tip:

    > How many SNVs should we be expecting?

    > Generally you can expect 1 SNV per 1000 bases and the exome of chromosome 22 is about 6-700kb (~2% of the 33500kb length of chr22). Therefore, we should expect to see ~700 SNVs.

    To further filter the dataset from above, we can use the *Filter and sort* tool.

    From the Galaxy tools panel, select
<div class=code>
Filter and Sort > Filter<br>
<br>
From the options:<br>
With following condition = c7>50 (filter on SNV quality (column 7)). <br>
This is a score generated by a combination of reads and base quality etc.<br>
<br>
keep other options as default and click execute<br>
</div>

    Note that the number of SNVs has gone down by ~95%, and we're down to 1059, which is not too far off the expected number.

    To visualize these variants we need to go back to IGV or any other genome browser of your choice.

    >  <img src="../media/tips.png" alt="Tip" height="42" width="42"/>
    >  Tip:
    > Interval file description<br>
    >* 1. Chromosome<br>
    >* 2. Start position (0-based)<br>
    >* 3. End position (1-based)<br>
    >* 4. Reference base at that position<br>
    >* 5. Coverage (# reads aligning over that position)<br>
    >* 6. Bases within reads<br>
    >* 7. Quality values (phred33 scale, see Galaxy wiki for more)<br>
    >* 8. Number of A variants<br>
    >* 9. Number of C variants<br>
    >* 10. Number of G variants<br>
    >* 11. Number of T variants<br>
    >* 12. Quality adjusted coverage<br>
    >* 13. Total number of deviants (if Convert coordinates to intervals? is set to yes)<br>

3.  Convert the variants file to [BED] format.
    * Start with renaming the variant file to e.g. 'NA12878.filtered.snps'
    * Change the filetype (under the pencil - see tips above) to bed.

    > Bed files are similar to Interval files, but follow a stricter structural format

    * Download the BED file by clicking on the disc icon and saving to local disk.
    * Assuming the alignment (bam file) is open, load the downloaded bed file.
    * Select *chr22* from the drop-down to see all SNVs.
    * The track - *NA12878.filtered.snps.bed* - now shows an zoomed out view of all SNVs.
    * Take a look again at the same region as earlier:
    * chr22:36,006,744-36,007,406

    <a href="../media/igv_snv.jpg"> <img src="../media/igv_snv.jpg" alt="IGV SNV view" width="640px" style="display:block; margin-left: auto; margin-right:auto;"/></a>
    > Is this a heterozygous variant?

    * Now zoom into the region chr22:35,947,503-35,947,667.
    > Is this a homozygous variant?

## 7. Calling small insertions and deletions (indels)

1.  Generate a pileup file from the aligned reads containing only indels.
    From the Galaxy tools panel, select

<div class=code>
NGS: SAM Tools > Generate pileup<br>
<br>
From the options:<br>
Select the BAM file to generate the pileup file for = sorted bam file<br>
Whether or not to print only output pileup lines containing indels = Print only lines containing indels<br>
Call consensus according to MAQ model? = yes<br>
<br>
Keep other options as default and click execute<br>
</div>

2.  Filtering pileup to get a list of indels

    Unlike the previous pileup, this file will only report putative indels. Like previously, we will need to filter the ~4000 indels down. Therefore, we will filter on:

    * removing poor quality bases (column 7)
    * the indel is covered by at least e.g. 10 reads (column 11)

    From the Galaxy tools panel, select
<div class=code>
Filter and Sort > Filter<br>
<br>
From the options:<br>
With following condition =  c7>50 and c11>20<br>
<br>
keep other options as default and click execute<br>
</div>

    This filtering strategy reduces the call set by 83% and we are left with **~700 small indels**.

    To visualize, we need to convert from *tabular* to *bed*. This is two step process.

    Firstly,
    > Click the pencil icon.

    > Under Attributes tab: make sure End column = 2

    > Under the Datatype tab: choose *Interval* and save

    Next, we can convert the *Interval* file to *BED* format.
    > Click the pencil icon.

    > Under Convert Format tab: choose *Convert Genomic Interval to BED*.

    > Rename this to indels.filtered


    Download this file and open it using IGV genome browser.

    Try looking at region chr22:31,854,409-31,854,460
    Can you see any indels?
    Also try looking at other regions by zooming out.  

<a href="../media/igv_indel.jpg"> <img src="../media/igv_indel.jpg" alt="IGV indel" width="640px" style="display:block; margin-left: auto; margin-right:auto;"/></a>

  ----------------------------------------------


[//]: # (These are reference links used in the body of this note and get stripped out when the markdown processor does it's job. There is no need to format nicely because it shouldn't be seen. Thanks SO - http://stackoverflow.com/questions/4823468/store-comments-in-markdown-syntax)

   [background]: <https://www.google.com/url?q=https://docs.google.com/document/pub?id%3D1NfythYcSrkQwldGMrbHRKLRORFFn-WnBm3gOMHwIgmE&sa=D&usg=AFQjCNE7C6wmK6Fiu-_ZJhc0RSBaxFSRbg>
   [1000 genomes]: <http://www.1000genomes.org/>
   [Melbourne Galaxy]: <http://galaxy-mel.genome.edu.au/galaxy>
   [FASTQ file]: <https://www.google.com/url?q=https://swift.rc.nectar.org.au:8888/v1/AUTH_a3929895f9e94089ad042c9900e1ee82/VariantDet_BASIC/NA12878.GAIIx.exome_chr22.1E6reads.76bp.fastq&sa=D&usg=AFQjCNFuof3Ud3BzcKkW54yL-1ySLIXPNg>
   [FASTQ]: <https://en.wikipedia.org/wiki/FASTQ_format>
   [BWA]: <http://bio-bwa.sourceforge.net/>
   [UCSC hg19]: <http://genome.ucsc.edu/cgi-bin/hgTracks>
   [SAM]: <https://samtools.github.io/hts-specs/SAMv1.pdf>
   [Pileup]: <https://www.google.com/url?q=https://docs.google.com/document/pub?id%3D1fouC29Lq0CXxQQCpuojrR5RXbdzMdxRf8ZID01XYNqI%23h.931a12a1a6ce&sa=D&usg=AFQjCNE-_ur_EqlMdu0A35ylXxlrNTlktA>
   [BED]: <https://genome.ucsc.edu/FAQ/FAQformat.html#format1>
