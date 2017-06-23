
<p>
<a href="https://www.melbournebioinformatics.org.au/"><img src="../../../img/melbioinf_logo.png" alt="Melbourne Bioinformatics logo" width="164"/></a>
<a href=http://genome.edu.au><img src="../media/gvl_logo.jpg" alt="GVL logo" align="right" width="112"/></a><br></br><br></br>
<br></br></p>
<p> </p>

# Variant Detection - Advanced Workshop

## Tutorial Overview

In this tutorial, we will look further at variant calling from sequence data. We will:

* Align NGS read data to a reference genome and perform variant calling, using somewhat different tools to those in the Basic workshop
* Carry out local realignment on our aligned reads
* Compare the performance of different variant calling tools
* Annotate our called variants with reference information

## Background

Some background reading and reference material can be found [here](var_detect_advanced_background.md).

The slides used in this workshop can be found [here](https://docs.google.com/presentation/d/1gsoS-ec1kCS7gBxLD4VAevhEMSsMA7TU7S8ZC_fzyEE/pub?start=false&loop=false&delayms=10000&slide=id.gd138ad700_0_0).

**Where is the data in this tutorial from?**

The data has been produced from human whole genomic DNA. Only reads that have mapped to a part of chromosome 20 have been used, to make the data suitable for an interactive tutorial. There are about one million 100bp reads in the dataset, produced on an Illumina HiSeq2000. This data was generated as part of the 1000 Genomes project: http://www.1000genomes.org/

## Preparation

1. **Make sure you have an instance of Galaxy ready to go.**
    * If you don't have your own - go to our [Galaxy-Tut](http://galaxy-tut.genome.edu.au/galaxy) or [Galaxy-Melbourne](http://galaxy-mel.genome.edu.au/galaxy) server.
    * Log in so that your work will be saved.
    If you don't already have an account on this server, select from the menu *User -> Register* and create one.

2. **Import data for the tutorial.**
    * We will import a pair of FASTQ files containing paired-end reads, and a
  VCF file of known human variants to use for variant evaluation.

    * **Method 1: Paste/Fetch data from a URL to Galaxy.**
    * In the Galaxy tools panel (left), under *BASIC TOOLS*, click on *Get Data* and choose *Upload File*.
    * Get the FASTQ files: click *Paste/Fetch data* and enter these URLs into the text box.
    If you put them in the same upload box, make sure there is a newline between the URLs so that they are really on separate lines.
<div class=code>
https://swift.rc.nectar.org.au:8888/v1/AUTH_377/public/variantCalling_ADVNCD/NA12878.hiseq.wgs_chr20_2mb.30xPE.fastq_1
<br>
https://swift.rc.nectar.org.au:8888/v1/AUTH_377/public/variantCalling_ADVNCD/NA12878.hiseq.wgs_chr20_2mb.30xPE.fastq_2
<br>
</div>
    * Select *Type* as **fastqsanger** and click *Start*. Note that you cannot use Auto-detect for the type here as there are different subtypes of FASTQ and Galaxy can't be sure which is which.
    * Get the VCF file: click *Paste/Fetch data* again to open a new text box, and paste the following URL into the box
<div class=code>
https://swift.rc.nectar.org.au:8888/v1/AUTH_377/public/variantCalling_ADVNCD/dbSNP135_excludingsitesafter129_chr20.vcf
<br>
</div>
    * This time, you can leave the *Type* on Auto-detect. Click *Start*.
    * Once the upload status for both sets of files turns *green*, you can click *Close*. You should now be able to see all three files in the Galaxy history panel (right).

    * **Method 2: Upload local data to Galaxy.** (In most cases, you won't need this for the tutorial)
    * *Use this method if you have your own files to upload, or if for any reason you find you need to manually download files for the tutorial.*
    * In the Galaxy tools panel (left), under *BASIC TOOLS*, click on *Get Data* and choose *Upload File*.
    * Click *Choose local file* and select the downloaded  FASTQ files. Select *Type* as **fastqsanger** and click *Start*.
    * Click *Choose local file* again and select the downloaded VCF file. Click *Start*.
    * Once the upload status for all files turns *green*, you can click *Close*. You should now be able to see all three files in the Galaxy history panel (right).

3. **Rename the datasets**
    * You should now have three files in your History, shown in the right-hand panel. If you used Method 1, the name of each dataset will be the full URL we got the file from. For convenience, we will give the datasets shorter names.
    * Click on the pencil icon to the top right of the dataset name (inside the green box) for the first dataset in your History. Note that the first dataset will be at the bottom! Shorten the name (you can just delete the first part) so that it is **NA12878.hiseq.wgs_chr20_2mb.30xPE.fastq_1**. Click *Save*.
    * Similarly, rename the second dataset to **NA12878.hiseq.wgs_chr20_2mb.30xPE.fastq_2**.
    * Similarly, rename the third dataset to **dbSNP135_excludingsitesafter129.chr20.vcf**.

## Section 1: Quality Control

The aim here is to evaluate the quality of the short data. If the quality is poor, then adjustments can be made - eg trimming the short reads, or adjusting your expectations of the final outcome!

1. **Analyse the quality of the reads in the FASTQ file.**
    * From the left hand tool panel in Galaxy, under *NGS ANALYSIS*, select *NGS: QC and manipulation -> FastQC*
    * Select one of the FASTQ files as input and *Execute* the tool.
    * When the tool finishes running, you should have an HTML file in your History. Click on the eye icon to view the various quality metrics.

Look at the generated FastQC metrics. This data looks pretty good - high per-base quality scores (most above 30).

## Section 2: Alignment and depth of coverage

In this step we map each of the individual reads in the sample FASTQ readsets to a reference genome, so that we will be able to identify the sequence changes with respect to the reference genome.

Some of the variant callers need extra information regarding the source of reads in order to identify the correct error profiles to use in their statistical variant detection model, so we add more information into the alignment step so that that generated BAM file contains the metadata the variant caller expects.

We will also examine the depth of coverage of the aligned reads across the genome, as a quality check on both the sequencing experiment and the alignment.

1. **Map/align the reads with Bowtie2 to the human reference genome.** We will use Bowtie2, which is one of several good alignment tools for DNA-seq data. Under *NGS ANALYSIS* in the tools panel, select the tool *NGS: Mapping -> Bowtie2*.
    * We have paired-end reads in two FASTQ files, so select *paired-end*.
    * Select the two FASTQ files as inputs.
    * Under *Select reference genome* select the human genome *hg19*.
    * Next we will add read group information. Read groups are usually used when we have reads from multiple experiments, libraries or samples, and want to put them into one aligned BAM file while remembering which read came from which group. In our case we only have one group, but the GATK tools need us to specify a read group in order to work correctly. Under *Set read groups information?* select *Set read groups (SAM/BAM specification)*. (Picard-style should also work.)
        * Set the read group identifier to "Tutorial_readgroup". This identifier needs to be a unique identifier for this read group. Since we only have one read group, it doesn't matter much what it is, but a common practice is to construct it out of information guaranteed to be unique, such as the library identifier plus Platform Unit (e.g. flowcell) identifier.
        * Set the sample name to "NA12878"
        * Set the platform to *ILLUMINA*
        * Set the library name to "Tutorial_library". Normally we would set this to identify the DNA library from our DNA extraction.
    * You can leave other read group information blank, and use default Bowtie2 settings. *Execute* the tool.
    * When the alignment has finished, you should rename the BAM file to something more convenient, such as **NA12878.chr20_2mb.30xPE.bam**.
    * *Note: we assume that you have seen BAM and SAM files before. If you have not you may want to try out the Basic Variant Calling workshop, or take the time now to convert your BAM file to a SAM file and examine the contents.*

2. **Visualise the aligned BAM file with IGV.** The Integrated Genome Viewer, IGV, is a very popular tool for visualising aligned NGS data. It will run on your computer (not on the server).
    * *Note: if you are already familiar with IGV, you may want to go through this section quickly, but it's still a good idea to launch IGV for use in later steps.*
    * In the green dataset box for your BAM file in the history panel, you will see some *display with IGV* links. Launch IGV by clicking the *web current* link. If IGV is already running on your computer, instead click the *local* link. If you have problems you can instead launch IGV by visiting https://www.broadinstitute.org/software/igv/download.
    * _If_ your BAM file was not automatically loaded, download and open it:
        * Download the BAM file AND the BAM index (BAI file) by clicking the floppy-disk icon in the green dataset window and selecting each file in turn. Make sure these two files are in the same directory.
        * In IGV, select the correct reference genome, *hg19*, in the top-left drop-down menu.
        * In IGV, open the BAM file using *File -> Load from File*.
    * Select chr20 in the IGV chromosomal region drop down box (top of IGV, on the left next to the organism drop down box).
    * Zoom in to the left hand end of chromosome 20 to see the read alignments - remember our reads only cover the first 2mb of the chromosome.
    * Scroll around and zoom in and out in the IGV genome viewer to get a feel for genomic data. Note that coverage is variable, with some regions getting almost no coverage (e.g. try chr20:1,870,686-1,880,895 - if you zoom right in to base resolution you’ll see that this region is very GC rich, meaning it’s hard to sequence. Unfortunately it also contains the first few exons of a gene...)

3. **Restrict the genomic region considered.** Later steps can be computationally intensive if performed on the entire genome. We will generate a genomic interval (BED) file that we will use to restrict further analyses to the first 2mb of chromosome 20, as we know our data comes from this region.
    * Under *BASIC TOOLS*, select the tool *Text manipulation -> Create single interval*. Enter these values:
        * Chromosome: chr20
        * Start position: 0
        * End position: 2000000
        * Name: chr20_2mb
        * Strand: plus
    * *Execute* this tool. This will create a small BED file specifying just one genomic region.
    * When the file is created, rename it to **chr20_2mb.bed**. Have a look at the contents of this BED file.

4. **Evaluate the depth of coverage of the aligned region.** Under *NGS COMMON TOOLSETS*, select the tool *NGS: GATK Tools 2.8 -> Depth of Coverage*.
    * Select the BAM file you just generated as the input BAM file.
    * Make sure the reference genome we aligned to is selected under *Using reference genome*.
    * Set *Output format* to *table*.
    * Restrict the analysis to only the region of interest: Set *Basic or Advanced GATK options* to *Advanced*. Click *Insert Operate on Genomic intervals* to add a new region and select the chr20_2mb.bed file from your history.
    * *Execute*.
    * This tool will produce a lot of files. We are most interested in the summaries. Examine the contents of:
        * **‘Depth of Coverage on data.... (output summary sample)’:** this file  will tell you the total depth of coverage for your sample across the genome (or in our case, across the first 2mb region we specified). It gives the total and mean coverage, plus some quantiles. This will give you an idea if there is something seriously wrong with your coverage distribution.
        * Your mean depth here should be ~24x. Note that ~89% of reference bases are covered by at least 15x coverage, which is a sort of informal agreed minimum for reasonable variant calling.
        * Also have a quick look at the **(per locus coverage)** file. It's not practical to go through this by hand, but you'll see that it gives coverage statistics for every site in the genome.
    * The other tables give you more detailed statistics on the level of coverage, broken down by regions etc. We don’t really need them so to keep our Galaxy history clean you can delete all the outputs of this step except for the ‘Depth of Coverage on data.... (output summary sample)’ file. Use the ‘X’ next to a history file to delete it.


## Section 3. Local realignment

Alignment to large genomes is a compromise between speed and accuracy. Since we usually have (at least) millions of reads, it becomes computationally too expensive to compare reads to one another - instead, high-throughput aligners such as Bowtie align each read individually to the reference genome.

This is often a problem where indels are present, as the aligner will be reluctant to align a read over an indel without sufficient evidence. This can lead to misalignment and to false positive SNVs. We can improve the performance of variant callers by carrying out a further step of ‘realignment’ after the initial alignment, considering all the reads in a particular genomic region collectively. In particular, this can provide enough collective evidence to realign reads correctly over suspected indels.

Both GATK and FreeBayes will benefit from local realignment around indels. Some tools, such as SamTools Mpileup, have alternative methods to deal with possible misalignment around indels.

If you skip this section, you can still carry out the variant calling steps in later sections by simply using the BAM file from Section 2.
However performing local realignment will improve the accuracy of our variant calls.

1. **Generate the list of indel regions to be realigned.** Under *NGS COMMON TOOLSETS*, select *NGS: GATK Tools 2.8 -> Realigner Target Creator*. This tool will look for regions containing potential indels.
    * Select your BAM file as input.
    * Select the correct reference genome (the genome used for alignment).
    * Restrict the analysis to only the region of interest: Set *Basic or Advanced GATK options* to *Advanced*. Click *Insert Operate on Genomic intervals* to add a new region and select the chr20_2mb.bed file from your history.
    * *Execute*.
    * If you examine the contents of the resulting intervals file, you will see a list of genomic regions to be considered in the next step.


2. **Realign the subsets of reads around the target indel areas.** Under *NGS COMMON TOOLSETS*, select *NGS: GATK Tools 2.8 -> Indel Realigner*. This step will realign reads and generate a new BAM file.
    * Select your BAM file as input.
    * Select the correct reference genome.
    * Under *Restrict realignment to provided intervals*, select the intervals file generated in the previous step.
    * *Execute*.
    * Rename the file to something easier to recognise, e.g. **NA12878.chr20_2mb.30xPE.realigned.bam**
    * To keep your history clean, you may want to delete the log files generated by GATK in the last two steps.


3. **Compare the realigned BAM file to original BAM around an indel.**
    * Open the realigned BAM file in IGV. You can use the *local* link to open it in your already-running IGV session and compare it to the pre-realignment BAM file.
    * Generally, the new BAM should appear identical to the old BAM except in the realigned regions.
    * Find some regions with indels that have been realigned (use the ‘Realigner Target Creator on data... (GATK intervals)’ file from the first step of realignment, it has a list of the realigned regions).
    * If you can’t find anything obvious, check region chr20:1163914-1163954;
    you should see that realignment has resulted in reads originally providing evidence of a ‘G/C’ variant at chr20:1163937 to be realigned with a 10bp insertion at chr20:1163835 and no evidence of the variant.


## Section 4. Calling variants with FreeBayes

FreeBayes is a Bayesian variant caller which assesses the likelihood of each possible genotype for each position in the reference genome, given the observed reads at that position, and reports back the list of variants (positions where a genotype different to homozygous reference was called). It  can be set to aggressively call all possible variants, leaving filtering to the user. FreeBayes generates a variant quality score (as do all variant callers) which can be used for filtering. FreeBayes will also give some phasing information, indicating when nearby variants appear to be on the same chromosome.

You can [read more about FreeBayes here](https://github.com/ekg/freebayes).

1. **Call variants with FreeBayes.** Under *NGS ANALYSIS*, select the tool *NGS: Variant Analysis -> FreeBayes*.
    * Select your realigned BAM file as input, and select the correct reference genome.
    * Under *Choose parameter selection level*, select "Simple diploid calling with filtering and coverage". This will consider only aligned reads with sufficient mapping quality and base quality. You can see the exact parameters this sets by scrolling down in the main Galaxy window to the Galaxy FreeBayes documentation section on *Galaxy-specific options*.
    * *Execute* FreeBayes.
    * When it has run, rename the resulting VCF file to something shorter, such as **NA12878.FreeBayes.chr20_2mb.vcf**.

2. **Check the generated list of variants**.
    * Click the eye icon to examine the VCF file contents.
    * How many variants exactly are in your list? (Hint: you can look at the number of lines in the dataset, listed in the green box in the History, but remember the top lines are header lines.)
    * What sort of quality scores do your variants have?
    * Open the VCF file in IGV using the dataset's *display in IGV local* link (using the *web current* link will open IGV again, and using *local* should use your already-running IGV). This will give an annotation track in IGV to visualise where variants have been called. Compare it to your BAM file.


## Section 5. Calling variants with GATK Unified Genotyper

For comparison, we will call variants with a second variant caller.

The GATK (genome analysis toolkit) is a set of tools from the Broad Institute. It includes the tools for local realignment, used in the previous step.

The GATK UnifiedGenotyper is a Bayesian variant caller and genotyper. You can also use the GATK HaplotypeCaller, which should be available on the GVL server you are using. It takes a little longer to run but is a more recent and sophisticated variant caller that takes into account human haplotype information. Both of these tools are intended primarily for calling diploid, germline variants.

You can [read more about the GATK here](http://www.broadinstitute.org/gatk).

1. **Call variants using Unified Genotyper.** Under *NGS COMMON TOOLSETS*, select the tool *NGS: GATK Tools 2.8 -> Unified Genotyper*.
    * Select your realigned BAM file as input, and select the correct reference genome.
    * UnifiedGenotyper can automatically label called variants that correspond to known human SNPs, if we provide it with a list of these. We have a VCF file of known SNPs which you imported as input data for this workshop. Under *dbSNP ROD file*, select the dataset *dbsnp135_excludingsitesafter129_chr20.vcf*.
    * Set *Genotype likelihoods calculation model to employ* to *SNP*.
    * Restrict the analysis to only the region of interest: Set *Basic or Advanced GATK options* to *Advanced*. Click *Insert Operate on Genomic intervals* to add a new region and select the chr20_2mb.bed file from your history.
    * *Execute*.
    * Rename the file to something useful eg **NA12878.GATK.chr20_2mb.vcf**.
    * The output file of interest is the VCF file. If you like, clean up your History by deleting the (log) and (metrics) files.


2. **Check the generated list of variants.**
    * Roughly how many variants are there in your VCF file (how many lines in the dataset?)
    * Click the eye icon to examine the contents of the VCF file. Notice that the *ID* column (the third column) is populated with known SNP IDs from the VCF file we provided.
    * Notice that the VCF header rows, and the corresponding information in the *INFO* column, is NOT the same as for the VCF file generated by FreeBayes. In general each variant caller gets to decide what information to put in the *INFO* column, so long as it adds header rows to describe the fields it uses.
    * Open the VCF file in IGV, using the link in the history panel.
    * Find a region where GATK has called a variant but FreeBayes hasn’t, and vice versa. Try chr20:1,123,714-1,128,378.


## Section 6. Evaluate variants

How can we evaluate our variants?

We've called variants on normal human DNA, so we expect to find variants with the typical characteristics of human germline variants. We know a lot about variation in humans from many empirical studies, including the 1000Genomes project, so we have some expectations on what we should see when we call variants in a new sample:

* We expect to see true variations at the rate of about 1 per 1000bp against the reference genome
* 85% of variations ‘rediscovered’ - that is, 85% already known and recorded in dbSNP (% dependent on the version of dbSNP)
* A transition/transversion (Ti/Tv) rate of >2 if the variants are high quality, even higher if the variants are in coding regions.

You can [read more about SNP call set properties here](http://genome.sph.umich.edu/wiki/SNP_Call_Set_Properties).

You may find that each of the variant callers has more variants than we would have expected - we would have expected around 2000 in our 2 megabase region but we see between 3000 and 5000. This is normal for variant calling, where most callers err on the side of sensitivity to reduce false negatives (missed SNVs), expecting the user to do further filtering to remove false positives (false SNVs).

We will also compare the output of our variant callers to one another - how many SNVs have they called in common? How many do they disagree on?

1. **Evaluate dbSNP concordance and Ti/Tv ratio using the GATK VariantEval tool.**
    * Under *NGS COMMON TOOLSETS*, select *NGS: GATK Tools 2.8 -> Eval Variants*.
    * Under *Input variant file*, select the VCF file you generated with FreeBayes.
    * Click *Insert Variant* to add a second input file, and under *Input variant file*, select the VCF file you generated with UnifiedGenotyper.
    * Set the reference genome to *hg19*.
    * Provide our list of known variants: make sure *Provide a dbSNP Reference-Ordered Data (ROD) file* is set to *Set dbSNP* and under *dbSNP ROD file* select the input reference variant file, **dbSNP135_excludingsitesafter129.chr20.vcf**.
    * We will avoid carrying out the full suite of comparisons for this tutorial, and look at a couple of metrics. Set *Basic or Advanced Analysis options* to *Advanced*. Then set the following:
        * *Eval modules to apply on the eval track(s)*:
            * *CompOverlap*
            * *TiTvVariantEvaluator*
        * *Do not use the standard eval modules by default*: check this option
    * *Execute*.


2. **Interpret the dbSNP concordance section of the evaluation report.**
    * Examine the contents of the ‘Eval Variants on data... (report)’
    * The first section of the report lists the overlap in variants between the generated VCF files and known human variation sites from dbSNP.
        * The *EvalRod* column specifies which of the input VCFs is analysed (input_0 = first VCF, input_1 = second etc). For us, these should be FreeBayes and GATK respectively.
        * The *CompRod* column specifies the set of variants against which the input VCFs are being compared; we have used dbsnp.
        * *Novelty*: whether the variants in this row have been found in the supplied dbSNP file or not (known = in dbSNP, novel = not in dbSNP). The rows containing *all* variants are the most informative summaries.
        * *nEvalVariants*: number of variant sites in EvalRod (i.e. all called variants).
        * *novelSites*: number of variant sites found in EvalRod but not CompRod (i.e. called variants not in dbSNP).
        * *CompOverlap*: number of variant sites found in both EvalRod and CompRod (i.e. called variants that ARE in dbSNP).
        * *compRate*: percentage of variant sites from EvalRod found in CompRod (=CompOverlap/nEvalVariants).
            * This metric is important, and it’s what people generally refer to as *dbSNP concordance*.
        * *nConcordant* and *concordantRate*: number and percentage of overlapping variants that have the same genotype.


3. **Interpret the TiTv section of the evaluation report.**
    * The second section of the report lists the transition/transversion ratio for different groups of variants from the callsets.
    * It generally follows the table format above. The most interesting metric is in the column labelled *tiTvRatio*.
    * The expectation is a Ti/Tv close to the TiTvRatioStandard of 2.38. Generally, the higher the better, as a high Ti/Tv ratio indicates that most variants are likely to reflect our expectations from what we know about human variation.
    * In this table, it can be useful to compare not just the rows labelled *all*, but also those labelled *known* and *novel*. Is the Ti/Tv ratio different for known dbSNP variants, as compared to novel variants in this individual?


4. **How much overlap is there in the call sets?**
    * One way to work out how many of the variants are in common between the call sets is to produce a Venn diagram to visualise the overlap. Since all our variants are on chromosome 20, we will do this by simply comparing the positions of the SNVs. So, we will first remove the VCF headers, leaving us just the variant rows, and then we will compare the variant coordinates.
    * Remove the header lines from a VCF file: select the tool *BASIC TOOLS -> Filter and Sort ->Select*.
        * As an input file, in *Select lines from*, select the VCF file you generated using FreeBayes.
        * Select *NOT Matching*.
        * As the pattern, enter "^#". *^* in regular expressions indicates the start of the line, so this is a regular expression that says we are detecting lines where there is a "#" character at the start of the line.
        * *Execute*. This should give you a new dataset, containing just the variant rows. Notice that the number of lines in this dataset now tells you how many variant calls you have!
    * Repeat the above steps to remove the header lines from VCF file that you generated using GATK UnifiedGenotyper.
    * Create a Venn diagram: select the tool *STATISTICS AND VISUALISATION -> Graph/Display Data -> proportional venn*.
        * Enter any title you like, e.g. "FreeBayes vs UnifiedGenotyper".
        * As *Input file 1*, select the first of the filtered files you just generated.
        * As *Column index*, enter 1. This is the second column, i.e. the column containing the position coordinate.
        * Under *as name*, enter a name for this input file, e.g. just "FreeBayes".
        * As *Input file 2*, select the second of the filtered VCF files.
        * Again as *Column index*, enter 1.
        * Under *as name*, enter a name for this input file, e.g. just "UnifiedGenotyper".
        * *Execute*.
    * You should get an HTML file containing a Venn diagram, containing the overlap of the two sets of SNV calls. The counts below show how many variants are in each part of the diagram:
        * "FreeBayes \ UnifiedGenotyper" indicates the difference, i.e. how many variants were called by FreeBayes and *not* UnifiedGenotyper.
        * "FreeBayes ∩ UnifiedGenotyper" indicates the intersection, i.e. how many variants were called by both variant callers.
    * You will probably find that FreeBayes calls variants more aggressively, and that there are more variants called by FreeBayes and not UnifiedGenotyper than vice versa. We could make FreeBayes calls more specific by filtering on e.g. variant quality score. We'd expect this to  remove many false positives, but also a few true positives.


## Section 7. Annotation

The variants we have detected can be annotated with information from known reference data. This can include, for instance

* whether the variant corresponds to a previously-observed variant in other samples, or across the population
* whether the variant is inside or near a known gene
* whether the variant is in a particular kind of genomic region (exon of a gene, UTR, etc)
* whether a variant is at a site predicted to cause a pathogenic effect on the gene when mutated

... and lots of other information!

Most human genes have multiple isoforms, i.e. multiple alternative splicings leading to alternative transcripts. The annotation information we get for a variant in the gene can depend on which transcript is used. In some cases, unless you request only one, you may see multiple alternative annotations for one variant.

For this workshop we will annotate our variants with the *SnpEff* tool, which has its own prebuilt annotation databases.

1. **Annotate the detected variants.** Select the tool *NGS ANALYSIS -> NGS: Annotation -> SnpEff*.
    * Choose any of your generated lists of variants as the input VCF file in the first field.
    * Select *VCF* as the output format as well. This will add the annotated information to the INFO column of the VCF file.
    * As *Genome source*, select *Named on demand*. Then as the genome, enter "hg19". This should work on any server, even if the SnpEff annotation database for hg19 has not been previously installed.
    * Make sure *Produce Summary Stats* is set to *Yes*.
    * *Execute* SnpEff.

2. **Examine the annotated information.** Click the eye icon on the resulting VCF file to view it.
    * You should see more information in the INFO column for each variant.
    * You should also see a few extra VCF header rows. In particular notice the new *INFO=<ID=EFF...* header row, listing the information added on the predicted effects of each variant.
    * Have a look through a few variants. Variants not inside or near known genes will not have much annotated information, and may simply have a short *EFF* field listing the variant as being in an "intergenic region". Variants in known functional regions of the genome will be annotated with much more information.
    * Try filtering to view only missense variants (i.e. substitutions which cause an amino acid change) by using the Galaxy tool *BASIC TOOLS -> Filter and Sort -> Select* on your annotated VCF file and filtering for lines matching the string "missense_variant".

3. **Examine the annotation summary stats.** The other file SnpEff produces is an HTML document with a summary of the annotations applied to the observed variants. Examine the contents of this file.
     * You will see summaries of types of variants, impact of variants, functional regions etc.
     Go through these tables and graphs and see if you understand what they represent.
     * Where are most of the variants found (in exons, introns etc)? Does this match what you'd expect?
