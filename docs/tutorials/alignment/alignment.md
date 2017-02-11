<style>


</style>
<div class="image-header-gvl">
    <a href=http://genome.edu.au><img src="../media/gvl_logo.jpg" alt="GVL logo" width="112"/></a>
    <a href="http://galaxyproject.org"><img src="../media/GalaxyLogoHighRes.png" width="164" /></a>
    <a href="https://www.melbournebioinformatics.org.au/"><img src="../../../img/melbioinf_logo.png" alt="Melbourne Bioinformatics logo" width="164"/></a>
</div>

# Alignment Tutorial

In this tutorial we will be performing some alignments of short reads to a longer reference (as outlined in earlier lectures.) It is split into two sections. You will perform the same analysis in both sections. The first is Alignment using the Galaxy bioinformatics workflow environment, the second is Alignment using the Unix/Linux command line. Both sections use the same tools.

We'll be using some data and steps from the Genomics Virtual Lab Variant Detection tutorials. You can follow these links if interested to see the full [Introductory Variant Detection](../variant_calling_galaxy_1/variant_calling_galaxy_1/) and [Advanced Variant Detection](../var_detect_advanced/var_detect_advanced/) tutorials.

## Section 1: Alignment using Galaxy

### Preparation

1. Open Galaxy using the IP address of your Cloud instance in Firefox or Chrome (whichever)
    * Unless you already have, create a log-in for yourself
    * On the top menu select: **User -> Register**
        * Enter your email address, choose a password, repeat it and add an (all lower case) one word name
        * Click Submit.
    * If you have already registered previously, just log in.

2. Create a new history and load the required data files
    * Click the <img src="../media/Galaxy-menu.png" width=20 /> icon in the History pane and **Create New**.

3. There's only one input file, so instead of importing a History, let's import the file directly:
    * From the tool panel, click on **Get Data -> Upload File**

        * Click on the **Paste/Fetch Data** button
          * Copy and paste the following web address into the URL/Text box:       *https://swift.rc.nectar.org.au:8888/v1/AUTH_377/public/variantCalling_BASIC/NA12878.GAIIx.exome_chr22.1E6reads.76bp.fastq*
          * Set the file format to *fastqsanger* (not fastqcsanger)
          * Click **Start**
    * Once the progress bar has reached 100%, click **Close**

4. Once the file is imported, you may want to give it a shorter name by clicking the <img src="../media/Galaxy-edit.png" width=20 /> icon on the dataset in the History pane.

5. Examine the input data.
    * In the History pane on the right, you should see a FASTQ file in the History. Click the <img src="../media/Galaxy-view.png" width=20 /> icon to view the contents.

The data for this workshop is the same as that used in the GVL Introductory Variant Detection tutorial. It is short read data from the exome of chromosome 22 of a single human individual. There are one million 76bp reads in the dataset, produced on an Illumina GAIIx from exome-enriched DNA. This data was generated as part of the 1000 Genomes project: [http://www.1000genomes.org/](http://www.1000genomes.org/).

### Alignment

#### Run BWA

For our alignment, we will use the tool BWA, which stands for "Burrows-Wheeler Aligner". You can see its website, if interested, at [http://bio-bwa.sourceforge.net/](http://bio-bwa.sourceforge.net/). We will be using the bwa mem variant.

Hopefully, everything will be working fine, so the first thing we need to do is run BWA.. To do this, use the following steps.

1. In the Galaxy tool pane on the left, under **NGS: Mapping** select the tool **Map with BWA-MEM**.

    * We want to align this data to the human genome. Have a look under **Using reference genome**. Choose *hg19*. This is Human reference genome 19 ([UCSC hg19](http://genome.ucsc.edu/cgi-bin/hgTracks?db=hg38&position=chr9%3A133252000-133280861&hgsid=469994169_9o1YboOfEP0gNexrwy41qevQrTdZ)).
    * Set **Single or Paired-end reads** to *Single*
    * Specify the FASTQ file as input - it may be automatically selected already.
    * For our alignment, it will be ok to use the default alignment parameters. However, have a look at the options you can change by setting **Select analysis mode** to *Full list of options*. Look through the parameters and see if you understand from lectures what they are for. You can also read documentation on these options by scrolling down the Galaxy BWA tool page.
    * Put the **Select analysis mode** back to *Simple Illumina Mode*
    * Click **Execute** and wait for your job to run.

This may take a few minutes and will produce a BAM file. While waiting you can have a look at the SAM/BAM format specification. Look at [this document](https://docs.google.com/document/pub?id=1fouC29Lq0CXxQQCpuojrR5RXbdzMdxRf8ZID01XYNqI#h.18e90b8fc68f), which gives you information on FASTQ, SAM, pileup and BED formats. You can also view the SAM/BAM format specification itself at [http://samtools.github.io/hts-specs/SAMv1.pdf](http://samtools.github.io/hts-specs/SAMv1.pdf).

Now we will convert the compressed BAM file to an uncompressed SAM file so we can see what it looks like.

1. In the Galaxy tool pane on the left, under **NGS: SAM tools** select the tool **BAM to SAM**.

    * For **BAM File to Convert** select your BAM file.
    * Click **Execute**.

Browse the resulting SAM file. To view the data click on the <img src="../media/Galaxy-view.png" width=20 /> icon next to the “BAM to SAM on … converted SAM” dataset in the History panel on the right-hand-side of the display. You should see header lines which begin with "@", followed by one alignment row per read. Refer to the docs above to understand what you're seeing.

One thing to notice is that the SAM file might be sorted in different ways, or might not be sorted at all. Two common ways of sorting a SAM or BAM file are:

* Alphabetically by the names of the reads. The read names are the first column of the SAM file, and come from the ">" lines in the FASTQ file.
* In order of genomic position, so that all the reads aligned to chromosome 1 are listed first (in coordinate order), followed by reads aligned to chromosome 2, etc. This information is in the third and fourth columns of the SAM file. This sort order is usually more useful for algorithms which need to access data that's aligned to a particular region of the genome.

## Section 2: View the alignment

### Load the VNC interface

We will use a Genome Browser called IGV to view the BAM alignment. It is installed on your cloud computer. To use it, we need to log in to your cloud computer's VNC interface.

1. Point your web browser at *http://your-ip-address*
2. This is the GVL dashboard. Half way down the page there is a link to the Lubuntu desktop. It is, *http://your-ip-address/vnc*, click on it.
3. You'll need to log in to the web page using the username: *Ubuntu*, password: *SummerCamp2016*. This will display the desktop of the cloud computer in your web browser. (You may need to reload the page for it to work properly.)
4. Login to the desktop as the *researcher* user (not ubuntu) using the same password.

Now we have logged into the cloud computer's visual desktop. You'll notice two icons on the desktop. A terminal and a shortcut to IGV.

### Download the data

The next thing we need to do is get the BAM files out of Galaxy and into the normal unix filesystem. To do this we need to load Galaxy within the VNC interface.

1. From the Linux menu button on the bottom left of the vnc interface, select **Internet->Firefox**
2. Point the resultant web browser (within the VNC interface) at *http://localhost/galaxy*
3. Login to this Galaxy interface via the **User** menu using the same login details you set up earlier. You should now see the history that you hve been working on.
4. Click on the file name of the BAM file in your history. This will expand it to show some additional details.
5. Click on the Download (<img src="../media/Galaxy-download.png" width=20 />) button. You'll need to download both the BAM file and BAM index file.

### View the BAM file in IGV

Now we will open IGV to view the bam file and browse the alignment.

1. Double click the IGV icon on the VNC desktop. Note that it will take some time to open.
2. We need to make sure we have the correct reference genome loaded that matches the one we used in our alignment. Make sure that the dropdown box on the top-left shows hg19.
3. Load the BAM you downloaded from Galaxy. Select **File -> Load from File** and supply the BAM file. (It is in the researcher user's Downloads directory.)

Now explore the visualisation of the BAM file!

Recall that our example data is from chromosome 22 only. In the contig drop-down menu, which says "All", instead select "22". Zoom in. If you're still having trouble finding the area the reads are mapped to, paste these coordinates into the box in the top bar of IGV: 22:17,432,499-17,502,283 . Then click "Go".

At the bottom, you should see an annotation track that is loaded into IGV by default, showing transcripts. Right click on the “Genes” track and select “Expanded”.

Zoom in to look at individual reads. Compare them to the gene annotation track, do you understand what is going on? Try zooming in far enough to see the reference genome sequence. Browse around and see if you can find some errors in the alignments. You can also explore IGV's features - try right-clicking on a read to change the display options.

## Section 3: Alignment using the command line

In this section of the tutorial we will perform the same alignment as before on the command line. We will use the same data and tools.

To use the command-line on your Research Cloud instance, ssh in to your instance as researcher@<your-ip-address> using ssh or Putty.

### Get the data

Once you are in the unix shell, you can get the input data for the tutorial using a command like `wget` or `curl`:

    wget https://swift.rc.nectar.org.au:8888/v1/AUTH_377/public/variantCalling_BASIC/NA12878.GAIIx.exome_chr22.1E6reads.76bp.fastq

This will download the datafile to the current directory.

### Set up the environment

Now you'll need to load the *bwa* and *samtools* modules. Try

    module avail

This command will show you a list of command line tool modules available to you. To use any of the tools, you need to load them into the current compute environment. You use the `module load` command for this. We need bwa and samtools so the commands are:

    module load bwa
    module load samtools/1.2

Then, run bwa with the command:	`bwa` to see some information on BWA usage. This will list BWA commands, such as *bwa index*, which indexes a reference genome, and *bwa mem*, which performs alignment. Type a command to see information on its usage. Try `bwa mem` to see how to run *bwa mem*.

Similarly, try just `samtools` to see a list of samtools commands.

### Perform the alignment

The Galaxy BWA tool actually wraps three commands; *bwa mem*, which performs alignment and produces a SAM file; *samtools view* which produces a BAM file from the SAM file and finally *samtools sort* which sorts the bam file so it is in reference position order.

We will replicate that process here.

How can you get an indexed reference genome? You could download a FASTA file and index it with `bwa index`. The human genome is fairly large and we won't have time to index it during the lab, but your instance should already have access to an indexed genome. Try

    ls galaxy_genomes
    ls galaxy_genomes/hg19
    ls galaxy_genomes/hg19/bwa_mem_index
    ls galaxy_genomes/hg19/bwa_mem_index/hg19/

Under galaxy_genomes are directories for different reference genomes, and each has pre-built indices for various tools including BWA. You may also need galaxy_genomes/hg19/sam_index when using some samtools commands.

You need to give BWA the prefix of the reference genome index files, which is *galaxy_genomes/hg19/bwa_index/hg19/hg19.fa*

So as the genome has already been indexed for us, we just need to point bwa at it.

So, now we can run bwa to align the reads to the reference.

We use the form of the command *bwa mem index_file reads_file > aligned_reads.sam*

To do this from your home directory, run the command:

    bwa mem galaxy_genomes/hg19/bwa_mem_index/hg19/hg19.fa NA12878.GAIIx.exome_chr22.1E6reads.76bp.fastq > aligned_reads.sam

This will load the reference index, align the reads and produce a SAM file called *aligned_reads.sam*.

Use `less` to examine it.

### Make the BAM file

Try to convert your SAM file to a BAM file using samtools view. The command you'll need is:

    samtools view -b -h aligned_reads.sam > aligned_reads.bam

You can look up the effect of all these command line flags by typing the view command by itself, or by reading the samtools documentation.

### Other SAM tools commands

Try some other samtools commands on your BAM file. Again, find out how to use them by typing them by themselves, or by reading the documentation. Some commands are:

  * *samtools sort* to sort the BAM file
  * *samtools index* to produce the index (.bai) file. You need this to visualise the BAM file in IGV.
  * *samtools flagstat* to get some basic information on mapped reads
  * *samtools view* to see the BAM files contents; this effectively converts it back to SAM

Now go back to the viewing BAM alignments section above. You'll need to sort the BAM file and create an index (.bai) file before you load the BAM file into IGV...

That's it, I hope you enjoyed the tutorial!
