<p>
<img src="../../../img/melbioinf_logo.png" alt="Melbourne Bioinformatics logo" align="left" width="40%"/><br />
<img src="../media/gvl_logo.jpg" alt="GVL logo" align="right" width="112"/>
</p>
<p><br></p>
<p><br></p>

# Introduction to Variant Calling using Galaxy

-----

## Tutorial Overview

In this tutorial we cover the concepts of detecting small variants (SNVs and indels) in human genomic DNA using a small set of reads from chromosome 22.

Note: The tutorial is designed to introduce the tools, datatypes and workflow of variation detection. We filter the variants manually to understand what is actually happening in variant calling. In practice the datasets would be much larger and you would carry out some extra steps to improve the quality of the called variants.

-----

## Learning Objectives

At the end of this tutorial you should:

 - Be familiar with the FASTQ format and base quality scores
 - Be able to align reads to generate a BAM file and subsequently generate a pileup file
 - Be able to run the FreeBayes variant caller to find SNVs and indels
 - Be able to visualise BAM files using the [Integrative Genomics Viewer (IGV)] and identify likely SNVs and indels by eye

-----

## Background

Some background reading material - [background]

#### Where is the data in this tutorial from?

The workshop is based on analysis of short read data from the exome of chromosome 22 of a single human individual. There are one million 76bp reads in the dataset, produced on an Illumina GAIIx from exome-enriched DNA. This data was generated as part of the [1000 genomes] Genomes project.

-----

## 1. Preparation

1. **Make sure you have an instance of Galaxy ready to go.**
    *  If you are not using your own Galaxy instance, you can use our [Galaxy Tutorial server] or [Galaxy Melbourne server].
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
        5. The dataset will have a very long name, as it's named after the full URL we got it from. Optionally, once the file is in your History, click the pencil icon in the upper-right corner of the green dataset box, then select the *Name* box and give the file a shorter name by removing the URL. Then click *Save*.

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

-----

## 2. Quality Control

The first step is to evaluate the quality of the raw sequence data. If the quality is poor, then adjustments can be made - e.g. trimming the short reads, or adjusting your expectations of the final outcome!

#### 1. Take a look at the FASTQ file

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

#### 2. Assessing read quality from the FASTQ files
1.  From the Galaxy tools panel, select `NGS: QC and manipulation > FastQC: Read Quality reports`.
<br>The input FASTQ file will be selected by default. Keep the other defaults and click execute.

   >  <img src="../media/tips.png" alt="Tip" height="42" width="42"/>
   Tip:
   Note the batch processing interface of Galaxy:  
   <ul>
   <li> grey = waiting in queue  
   <li> yellow = running  
   <li> green = finished  
   <li> red = tried to run and failed
   </ul>

When the job has finished, click on the eye icon to view the newly generated data (in this case a set of quality metrics for the FASTQ data).
Look at the various quality scores. The data looks pretty good - *high Per base sequence quality* (avg. above 30).

*Note that the 'Sequence Duplication Levels' are marked as high. Normally we would run another tool to remove duplicates (technical PCR artifacts) but for the sake of brevity we will omit this step.*

-----

## 3. Alignment to the reference - (FASTQ to BAM)

The basic process here to map individual reads - from the input sample FASTQ file - to a matching region on the reference genome.

#### 1. Align the reads with BWA
1.  Map/align the reads with the [BWA] tool to Human reference genome 19 (hg19) [UCSC hg19].
    From the Galaxy tools panel, select
<div class=code>
*NGS: Mapping > Map with BWA-MEM [3-5mins]*<br>
<br>
From the options:<br>
Using reference genome: set to hg19<br>
Single or Paired-end reads: set to Single<br>
<br>
Make sure your fastq file is the input file.
<br>
Keep other options as default and click execute.<br>
</div>

    *Note: This is the longest step in the workshop and will take a few minutes, possibly more depending on how many people are also scheduling mappings*

2.  Sort the BAM file.
  From the Galaxy tools panel, select
<div class=code>
*NGS: SAM Tools > Sort BAM dataset*
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

#### 2. Examine the alignment
1.  To examine the output sorted BAM file, we need to first convert it into readable [SAM] format.
  From the Galaxy tools panel, select
<div class=code>
*NGS: SAM Tools > BAM-to-SAM*
<br><br>
From the options:
<br>
BAM File to Convert: set to the output to the sorted BAM file
<br><br>
Keep other options as default and click execute
<br>
</div>

2.  Examine the generated Sequence Alignment Map (SAM) file.
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


#### 3. Assess the alignment data

We can generate some mapping statistics from the BAM file to assess the quality of our alignment.

1.  Run IdxStats. From the Galaxy tools panel, select
    <div class=code>
    *NGS: SAM Tools > IdxStats*
    <br>
    Select the sorted BAM file as input. Keep other options as default and click execute.
    <br>
    </div>

    IdxStats generates a tab-delimited output with four columns.
    Each line consists of a reference sequence name (e.g. a chromosome),
    reference sequence length, number of mapped reads and number of placed but
    unmapped reads.

    We can see that most of the reads are aligning to chromosome 22 as expected.

2.  Run Flagstat. From the Galaxy tools panel, select
    <div class=code>
    *NGS: Sam Tools > Flagstat*
    <br><br>
    From the options:
    <br>
    The BAM: select the sorted BAM file
    <br><br>
    Keep other options as default and click execute
    <br>
    </div>

    Note that in this case the statistics are not very informative. This is because the dataset has been generated for this workshop and much of the noise has been removed (and in fact we just removed a lot more noise in the previous step); also we are using single ended read data rather than paired-end so some of the metrics are not relevant.

-----

## 4. Visualise the BAM file.

To visualise the alignment data:

1.  Click on the sorted BAM file dataset in the History panel.
2.  Click on "Display with IGV **web current**". This should download a .jnlp
    Java Web Start file to your computer. Open this file to run IGV.
    (You will need Java installed on your computer to run IGV). NOTE: If IGV is already open on your computer, you can click "**local**" instead of "web current", and this will open the BAM file in your current IGV session.
3.  Once IGV opens, it will show you the BAM file. This may take a bit of time as the data is downloaded.
4.  Our reads for this tutorial are from chromosome 22, so select **chr22** from the second
    drop box under the toolbar. Zoom in to view alignments of reads to the reference genome.

<a href="../media/igv1.jpg"><img src="../media/igv1.jpg" alt="IGV view 1" width="640px" style="display:block; margin-left: auto; margin-right:auto;"/></a>

Try looking at region `chr22:36,006,744-36,007,406`

Can you see a few variants?  

<a href="../media/igv_mb.jpg"><img src="../media/igv_mb.jpg" alt="IGV view 1" width="640px" style="display:block; margin-left: auto; margin-right:auto;"/></a>

Don't close IGV yet as we'll be using it later.

-----

## 5. Generate a pileup file

A pileup is essentially a column-wise representation of the aligned reads - at the base level - to the reference. The pileup file summarises all data from the reads at each genomic region that is covered by at least one read. Each row of the pileup file gives similar information to a single vertical column of reads in the IGV view.

The current generation of variant calling tools do not output pileup files, and you don't need to do this section in order to use FreeBayes in the next section. However, a pileup file is a good illustration of the evidence the variant caller is looking at internally, and we will produce one to see this evidence.

1. **Generate a pileup file:**

    From the Galaxy tools panel, select
    <div class=code>
    *NGS: SAMtools > Generate Pileup*<br>
    <br>
    From the options:<br>
    Call consensus according to MAQ model = Yes<br>
    This generates a called 'consensus base' for each chromosomal position.<br>
    <br>
    Keep other options as default and click execute<br>
    </div>

    For each output file, Galaxy tries to assign a datatype attribute to every file. In this case, you'll need to manually assign the datatype to *pileup*.

    1. First, rename the output to something more meaningful by clicking on the pencil icon
    2. To change the datatype, click on the Datatype link from the top tab while you're editing attributes.
    3. For downstream processing we want to tell Galaxy that this is a *pileup* file. From the drop-down, select Pileup and click save.

    >  <img src="../media/tips.png" alt="Tip" height="42" width="42"/>
    >  Tip:
    >The pileup file we generated has 10 columns:<br>
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

    > Further information on (9):<br>
    > Each character represents one of the following (the longer this string, higher the coverage):
    >   <ul>
    >   <li> . = match on forward strand for that base
    >   <li> , = match on reverse strand
    >   <li> ACGTN = mismatch on forward
    >   <li> acgtn = mismatch on reverse
    >   <li> +[0-9]+[ACGTNacgtn]+' = insertion between this reference position and the next
    >   <li> -[0-9]+[ACGTNacgtn]+' = deletion between this reference position and the next
    >   <li> ^ = start of read
    >   <li> $ = end of read
    >   <li> BaseQualities = one character per base in ReadBases, ASCII encoded Phred scores
    >   </ul>

2. **Filter the pileup file:**

    If you click the eye icon to view the contents of your pileup file, you'll see that the visible rows of the file aren't very interesting as they are outside chromosome 22 and have very low coverage. Let's filter to regions with coverage of at least 10 reads.
      <br>
      From the Galaxy tools panel, select:
      <div class=code>
      *NGS: SAM Tools > Filter Pileup*<br>
      <br>
      From the options:<br>
      which contains = Pileup with ten columns (with consensus)<br>
      Do not report positions with coverage lower than = 10<br>

    Try this filtering step two ways:

      * First, set the parameter *Only report variants?* to **No**. This will give you all locations with a coverage of at least 10 and sufficient read quality. This is similar to the information you see when you look at IVG: each pilup row corresponds to a column of aligned bases at one genomic location.
      * Then, repeat the step but set *Only report variants?* to **Yes**. This will effectively do variant calling: it will give you only locations that have some evidence that there might be a variant present. This variant calling is not very stringent, so you will still get lots of rows. We could filter further to, for instance, variants with high quality scores.

    Examine your two pileup files and understand the difference between them. Which coordinates are present in each? What do the bases look like in one compared to the other? Compare the variant quality score (in column 6) to the bases listed on each row.

## 6. Call variants with FreeBayes

FreeBayes is a Bayesian variant caller which assesses the likelihood of each possible genotype for each position in the reference genome, given the observed reads at that position, and reports back the list of possible variants. We look at it in more detail in the [Advanced Variant Calling](../var_detect_advanced/var_detect_advanced) tutorial.

1. **Call variants with FreeBayes.** Under *NGS ANALYSIS*, select the tool *NGS: Variant Analysis -> FreeBayes*.
    * Select your sorted BAM file as input, and select the correct reference genome.
    * Under *Choose parameter selection level*, select "Simple diploid calling with filtering and coverage". This will consider only aligned reads with sufficient mapping quality and base quality. You can see the exact parameters this sets by scrolling down in the main Galaxy window to the Galaxy FreeBayes documentation section on *Galaxy-specific options*.
    * *Execute* FreeBayes.

2. **Check the generated list of variants**.
    * Click the eye icon to examine the VCF file contents. The VCF format is described below - make sure you can identify the header rows and the data, and understand the important columns.
    * How many variants exactly are in your list? (Hint: you can look at the number of lines in the dataset, listed in the green box in the History, but remember the top lines are header lines.)
    * What sort of quality scores do your variants have?

    FreeBayes, like most variant callers, produces a [Variant Call Format (VCF)](http://vcftools.sourceforge.net/specs.html) file.

    VCF consists of a header section and a data section. The header section has some information about the file and the parameters used to produce it. The header also specifies what information is stored in the INFO, FILTER and FORMAT columns, as this is different for different variant callers.

    The data section has several columns. For this tutorial, you should concentrate on CHROM, POS, REF, ALT and QUAL. CHROM and POS describe the variant's genomic location, REF and ALT describe the variant's nucleotides relative to the reference, and QUAL is a quality score giving FreeBayes' confidence in the correctness of this variant call.

    The columns in more detail are:

    Col      | Field     | Description
    ---------|-----------|-------------------------------------------------------------------------------------
    1        |CHROM      |Chromosome name
    2        |POS        |1-based position. For an indel, this is the position preceding the indel.
    3        |ID         |Variant identifier (optional). Usually the dbSNP rsID.
    4        |REF        |Reference sequence at POS involved in the variant. For a SNP, it is a single base.
    5        |ALT        |Comma delimited list of alternative sequence(s) seen in our reads.
    6        |QUAL       |Phred-scaled probability of all samples being homozygous reference.
    7        |FILTER     |Semicolon delimited list of filters that the variant fails to pass.
    8        |INFO       |Semicolon delimited list of variant information.
    9        |FORMAT     |Colon delimited list of the format of individual genotypes in the following fields.
    10+      |Sample(s)  |Individual genotype information defined by FORMAT.

    For even more detail on VCF files, you can look at the [*VCF format specification*](http://vcftools.sourceforge.net/specs.html).

3. **Visualise the variants and compare files**
    * Open the VCF file in IGV using the dataset's *display in IGV local* link (using the *web current* link will open IGV again, and using *local* should use your already-running IGV). This will give an annotation track in IGV to visualise where variants have been called. Compare it to your BAM file.
    * Take a look again at the same region as earlier: `chr22:36,006,744-36,007,406`
    * Try comparing to the corresponding location in the pileup file. You can filter to the same window as we just opened in IGV with the tool *Filter and Sort > Filter*. Choose your previously-filtered pileup file as input, and set the filter condition to `c1=="chr22" and c2 > 36006744 and c2 < 36007406`.

4. **Optional: filter variants**: See if you can work out how to filter your VCF file to variants with quality scores greater than 50. You can use the *Filter and Sort: Filter* tool we used above.

## 7. Further steps

We've seen how to:

 - Align the raw data (sequence reads) to a reference genome
 - Generate variant calls from aligned reads
 - Interpret the various file formats used in storing reads, alignments and variant calls
 - Visualise the data using IGV

For real variant calling, you will probably want to carry out clean-up steps on your BAM file to improve the quality of the calls, and do further filtering and selection on the resulting variants.

We look at some further steps in the [Advanced Variant Calling](../var_detect_advanced/var_detect_advanced) tutorial.


  ----------------------------------------------


[//]: # (These are reference links used in the body of this note and get stripped out when the markdown processor does its job. There is no need to format nicely because it shouldn't be seen. Thanks SO - http://stackoverflow.com/questions/4823468/store-comments-in-markdown-syntax)

   [background]: <https://www.google.com/url?q=https://docs.google.com/document/pub?id%3D1NfythYcSrkQwldGMrbHRKLRORFFn-WnBm3gOMHwIgmE&sa=D&usg=AFQjCNE7C6wmK6Fiu-_ZJhc0RSBaxFSRbg>
   [1000 genomes]: <http://www.1000genomes.org/>
   [Galaxy Tutorial server]: <http://galaxy-tut.genome.edu.au/galaxy>
   [Galaxy Melbourne server]: <http://galaxy-mel.genome.edu.au/galaxy>
   [FASTQ file]: <https://www.google.com/url?q=https://swift.rc.nectar.org.au:8888/v1/AUTH_a3929895f9e94089ad042c9900e1ee82/VariantDet_BASIC/NA12878.GAIIx.exome_chr22.1E6reads.76bp.fastq&sa=D&usg=AFQjCNFuof3Ud3BzcKkW54yL-1ySLIXPNg>
   [FASTQ]: <https://en.wikipedia.org/wiki/FASTQ_format>
   [BWA]: <http://bio-bwa.sourceforge.net/>
   [UCSC hg19]: <http://genome.ucsc.edu/cgi-bin/hgTracks>
   [SAM]: <https://samtools.github.io/hts-specs/SAMv1.pdf>
   [Pileup]: <https://www.google.com/url?q=https://docs.google.com/document/pub?id%3D1fouC29Lq0CXxQQCpuojrR5RXbdzMdxRf8ZID01XYNqI%23h.931a12a1a6ce&sa=D&usg=AFQjCNE-_ur_EqlMdu0A35ylXxlrNTlktA>
   [BED]: <https://genome.ucsc.edu/FAQ/FAQformat.html#format1>
   [Integrative Genomics Viewer (IGV)]: <http://software.broadinstitute.org/software/igv/>
