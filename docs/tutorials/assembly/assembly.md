
<p>
<a href=http://vlsci.org.au><img src="media/vlsci_logo.jpg" alt="VLSCI logo" align="left" width="164"/></a>
<a href=http://genome.edu.au><img src="media/gvl_logo.jpg" alt="GVL logo" align="right" width="112"/></a>
</p>
<p></p>

# Microbial de novo Assembly for Illumina Data

## *Introductory Tutorial*

<img src="media/D4D_5503X.jpg" width="70"/>

Written and maintained by [Simon Gladman](mailto:simon.gladman@unimelb.edu.au) - VLSCI

<!-- toc -->
## Contents

1. [Tutorial Overview](#1-tutorial-overview)
2. [Background](#2-background-15-min) [15 min]
3. [Preparation](#3-preparation-15-min) [15 min]
4. [Section 1: Quality control](#4-section-1-quality-control-30-mins) [30 mins]
5. [Section 2: Assemble reads into contigs with Velvet and the Velvet Optimiser](#5-section-2-assemble-reads-into-contigs-with-velvet-and-the-velvet-optimiser-45-min) [45 min]
6. [Section 3. Extension](#6-section-3-extension-20-min) [20 min]
7. [References](#7-references)

## Tutorial Overview

In this tutorial we cover the concepts of Microbial de novo assembly using a very small synthetic dataset from a well studied organism.


**What’s not covered**

This tutorial covers the basic aspects of microbial de novo assembly from Illumina paired end or single end reads.
It does not cover more complicated aspects of assembly such as:
* Incorporation of other raw data types (454 reads, Sanger reads)
* Gap filling techniques for “finishing” an assembly
* Measuring the accuracy of assemblies

## Background [15 min]

Read the [background to the workshop here](assembly_background.md)

**Where is the data in this tutorial from?**

The data for this tutorial is from a whole genome sequencing experiment of a multi-drug resistant strain of the bacterium Staphylococcus aureus. The DNA was sequenced using an Illumina GAII sequencing machine. The data we are going to use consists of about 4 million x 75 base-pair, paired end reads (two FASTQ read files, one for each end of a DNA fragment.) The data was downloaded from the [NCBI Short Read Archive (SRA)](http://www.ncbi.nlm.nih.gov/sra/) (http://www.ncbi.nlm.nih.gov/sra/). The specific sample is a public dataset published in April 2012 with SRA accession number ERR048396.

We will also use a FASTA file containing the sequences of the Illumina adapters used in the sequencing process. It is desirable to remove these as they are artificial sequences and not part of the bacterium that was sequenced.

We will use software called [Velvet](https://www.ebi.ac.uk/~zerbino/velvet/) (Zerbino et al 2008) for the main de novo assembly, as well as some other peripheral software for pre- and post-processing of the data. Details of these can be found in the background document linked above.

**The protocol:**

We are performing a de novo assembly of the read data into contigs and then into scaffolds (appropriately positioned contigs loosely linked together). We firstly need to check the quality of the input data as this will help us choose the most appropriate range of input parameters for the assembly and will guide us on an appropriate quality trimming/cleanup strategy. We will then use an iterative method to assemble the reads using the [Velvet Optimiser](http://www.vicbioinformatics.com/software.velvetoptimiser.shtml) (a program that performs lots of Velvet assemblies searching for an optimum outcome.) Once this is complete we will obtain summary statistics on the final results (contigs) of the assembly.

Follow this [link for an overview of the protocol](protocol.md)

**The protocol in a nutshell:**

*Input:* Raw reads from sequencer run on microbial DNA sample.

*Output:* File of assembled scaffolds/contigs and associated information.


## Preparation [15 min]

### Login to Galaxy

1. Open a browser and go to a Galaxy server. (what is [Galaxy](assembly_background.md)?)
  * You can use a galaxy server of your own or
  * [Galaxy Tute](http://galaxy-tut.genome.edu.au) at genome.edu.au

2. Register as a new user if you don’t already have an account on that particular server



  >  <img src="media/tips.png" alt="Tip" height="42" width="42"/>
  NOTE: Firefox/Safari/Chrome all work well, Internet Explorer not so well.

### Import the DNA read data for the tutorial.

**You can do this in a few ways. If you're using [galaxy-tut.genome.edu.au](http://galaxy-tut.genome.edu.au):**
  1. Go to **Shared Data -> Published Histories** and click on *‘Microbial_assembly_input_data’*. Then click **'Import History'** at top right, wait for the history to be imported to your account, and then **‘start using this history’**.
  2. This will create a new Galaxy history in your account with all of the required data files.
  3. Proceed to step 4.

**If you are using a different Galaxy server, you can upload the data directly to Galaxy using the file URLs.**
  1. On the Galaxy tools panel, click on **Get data -> Upload File**.
  2. Click on the **Paste/Fetch Data** button.
  3. Paste the URL: https://swift.rc.nectar.org.au:8888/v1/AUTH_377/public/Assembly/ERR048396_1.fastq.gz into the text box. Change the type to *fastqsanger* (Not *fastqcsanger*).
  4. Click on the **Paste/Fetch Data** button again.
  5. Paste the URL: https://swift.rc.nectar.org.au:8888/v1/AUTH_377/public/Assembly/ERR048396_2.fastq.gz into the text box and change it's type to *fastqsanger* as well.
  6. Repeat the process for the last URL: https://swift.rc.nectar.org.au:8888/v1/AUTH_377/public/Assembly/illumina_adapters.fna, but make it's type *fasta*
  7. Click on the **Start** button. Once all of the uploads are at 100%, click on the **Close** button.
  8. When the files have finished uploading, rename them to ‘ERR048396_1.fastq’, ‘ERR048396_2.fastq’ and ‘illumina_adapters.fna’ respectively by clicking on the <img src="media/Galaxy-edit.png" width=20 /> icon to the top right of the file name in the right hand Galaxy panel (the history panel)


You should now have the following files in your Galaxy history:

* *ERR048396_1.fastq* - forward reads in fastq format
* *ERR048396_2.fastq* - reverse reads in fastq format
* *illumina_adapters.fa* - Illumina adapter sequences in fasta format


### View the fastq files
Click on the <img src="media/Galaxy-view.png" width=20 />  icon to the top right of each fastq file to view the first part of the file
If you’re not familiar with the FASTQ format, click here for an overview

NOTE: If you log out of Galaxy and log back in at a later time your data and results from previous experiments will be available in the right panel of your screen called the ‘History.’

----------------------------------------------

## Section 1: Quality control [15 mins]

The basic process here is to collect statistics about the quality of the reads in the sample FASTQ readsets. We will then evaluate their quality and choose an appropriate regime for quality filtering using Trimmomatic (a fastq read quality trimmer.)

[More detailed description of FastQC quality analysis can be found here.](https://docs.google.com/document/pub?id=16GwPmwYW7o_r-ZUgCu8-oSBBY1gC97TfTTinGDk98Ws)

[More detailed description of Trimmomatic read quality filtering can be found here.](http://www.usadellab.org/cms/index.php?page=trimmomatic)

### Run FastQC on both input read files
1. From the tools menu in the left hand panel of Galaxy, select **NGS QC and manipulation > FastQC: Comprehensive QC** (down the bottom of this category) and run with these parameters:
  * FASTQ reads: *ERR048396_1.fastq*
  * Use default for other fields
  * Click **Execute**

Screenshot of this process can be seen [here.](assembly-background.md)

> Note: This may take a few minutes, depending on how busy Galaxy is.

2. Now repeat the above process on the second read file: *ERR048396_2.fastq*

It is important to do both read files as the quality can be very different between them.
Examine the FastQC output
You should have two output objects from the first step:
FastQC_ERR048396_1.fastqc.html
FastQC_ERR048396_2.fastqc.html
These are a html outputs which show the results of all of the tests FastQC performed on the read files.
Click on the “eye” of each of these objects in turn to see the FastQC output.
The main parts of the output to evaluate are:
Basic statistics. This section tells us that the ASCII quality encoding format used was Sanger/Illumina 1.9 and the reads are length 75 and the percent GC content of the entire file is 35%.
Per base sequence quality. In the plot you should see that most of the early bases are up around the '32' mark and then increase to 38-40, which is very high quality; The spread of quality values for the last few bases increases and some of the outliers have quality scores of less than 30. This is a very good quality dataset. 20 is often used as a cutoff for reliable quality.
Screenshot can be seen here
Quality trim the reads using Trimmomatic.
From the tools menu in the left hand panel of Galaxy, select NGS QC and manipulation > Trimmomatic (down the bottom of this category) and run with these parameters (only the non-default selections are listed here):
Direction 1 fastq reads to trim: ERR048396_1.fastq
Direction 1 fastq reads to trim: ERR048396_2.fastq
Quality encoding: phred33
Clip Illumina adapters? checkbox: checked
Fasta of adapters to clip: illumina_adapters.fa
Trim leading bases, minimum quality: 15
Trim trailing bases, minimum quality: 15
Minimum length read: 35
Execute
Screenshot of this process can be seen here.
Examine the Trimmomatic output FastQ files.
You should have 3 objects from the output of Trimmomatic:
Trimmomatic on data 3, data 1, and data 2: Dir1 trimmed pairs
Trimmomatic on data 3, data 1, and data 2: Dir2 trimmed pairs
Trimmomatic on data 3, data 1, and data 2: trimmed reads
Click on the “eye” symbol on one of the objects to look at its contents. You’ll notice that not all of the reads are the same length now, as they have had the illumina adapters cut out of them and they’ve been quality trimmed.
Completed Galaxy history for this section (in SharedData>Published Histories): Microbial_assembly_section1


----------------------------------------------

[//]: # (These are reference links used in the body of this note and get stripped out when the markdown processor does it's job. There is no need to format nicely because it shouldn't be seen. Thanks SO - http://stackoverflow.com/questions/4823468/store-comments-in-markdown-syntax)
