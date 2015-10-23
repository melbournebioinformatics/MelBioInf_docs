
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
##Contents

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

**Login to Galaxy**

1. Open a browser and go to a Galaxy server. (what is [Galaxy](assembly_background.md)?)
  * You can use a galaxy server of your own or
  * [Galaxy Tute](http://galaxy-tut.genome.edu.au) at genome.edu.au

  >  <img src="media/tips.png" alt="Tip" height="42" width="42"/>
  Tip:
  NOTE: Firefox/Safari/Chrome all work well, Internet Explorer not so well
  Register as a new user: User>Register or login if you already have an account

2. Register as a new user if you don’t already have an account on that particular server


Import the DNA read data for the workshop. You can do this in a few ways, of which by far the easiest is:
Go to shared Data -> Published Histories and click on  ‘Microbial_assembly_input_data’. Then click 'Import History' at top right, wait for the history to be imported to your account, and then ‘start using this history’.
This will create a new Galaxy history in your account with all of the required data files
Proceed to step 3.

Alternately:
You can download the data from the GVL repository using the following URLs, then upload it to Galaxy for use.
Download the following files:
https://swift.rc.nectar.org.au:8888/v1/AUTH_377/public/Assembly/ERR048396_1.fastq.gz
https://swift.rc.nectar.org.au:8888/v1/AUTH_377/public/Assembly/ERR048396_2.fastq.gz
https://swift.rc.nectar.org.au:8888/v1/AUTH_377/public/Assembly/illumina_adapters.fna
Once files are downloaded, for each of them:
In Galaxy tools panel click on Get data>Upload File
Click 'Choose File', locate the local copy file and upload
For the ‘.fastq.gz’ files, make sure you select ‘fastqsanger’ under the file format to tell Galaxy this is a fastq file. Galaxy will unzip it automatically. For the ‘.fna’ file, make sure you select ‘fasta’.
Click ‘Execute’ to upload. (This may take a couple of minutes for each file.)
Or:
You can upload the data directly to Galaxy using the file URLs.
In Galaxy tools panel click on Get data>Upload File
copy the URLs from the previous step into the URL/Text box one at a time, and Galaxy will download and import them directly. Remember to select the ‘fastqsanger’ File Format for the sequence files and ‘fasta’ for the ‘.fna’ file.
When the files have finished uploading, rename them to ‘ERR048396_1.fastq’, ‘ERR048396_2.fastq’ and ‘illumina_adapters.fna’ respectively by clicking on the pencil icon to the top right of the file name in the right hand Galaxy panel (the history panel)


You should now have the following files in your Galaxy history:
A set of paired end reads consisting of 2 files:
ERR048396_1.fastq
ERR048396_2.fastq
Illumina adapter list as illumina_adapters.fa


View the fastq files
Click on the eye icon to the top right of each fastq file to view the first part of the file
If you’re not familiar with the FASTQ format, click here for an overview

NOTE: If you log out of Galaxy and log back at a later time your data and results from previous experiments will be available in the right panel of your screen called the ‘History’


Completed Galaxy history for this section (in SharedData>Published Histories): Microbial_assembly_complete



  ----------------------------------------------

[//]: # (These are reference links used in the body of this note and get stripped out when the markdown processor does it's job. There is no need to format nicely because it shouldn't be seen. Thanks SO - http://stackoverflow.com/questions/4823468/store-comments-in-markdown-syntax)
