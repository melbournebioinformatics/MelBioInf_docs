<img src="media/melbioinf_logo.png" width="350"> <img src="media/PRIMARY_A_Vertical_Housed_RGB.png" width="150">

# Introduction to Genome Browsers

Anticipated workshop duration when delivered to a group of participants is **2 hours**.  

For queries relating to this workshop, contact Melbourne Bioinformatics (bioinformatics-training@unimelb.edu.au).

## Overview
This tutorial will introduce you to the genome browsers format and illustrate how some freely avalable genome browsers can be used to interrogate a variety of data types.

### Topic

* [x] Genomics
* [x] Transcriptomics
* [ ] Proteomics
* [ ] Metabolomics
* [ ] Statistics and visualisation
* [ ] Structural Modelling
* [x] Basic skills


### Skill level

* [x] Beginner  
* [ ] Intermediate  
* [ ] Advanced  

This workshop is designed for participants with no previous experience of using Genome Browsers and no programming experience.


### Description

*learn how to make the most of Genome Browsers !*

Genome browsers are invaluable for viewing and interpreting any data that can be anchored to a position on a genome. These include gene location, genomic variation, transcription, and the many types of regulatory data such as methylation sites and transcription factor binding sites. Genome browsers enabling the viewing of one type of data in the context of other data types which can reveal important information about mutations and gene expression in normal development and disease.

This tutorial is in two parts.</br>
Section 1 introduces general features of genome browsers. </br>
Section 2 has a number of hands on tutorials to get you used to using different genome browsers:
* [UCSC genome Browser](https://genome.ucsc.edu/)
* [Integrative Genomics Viewer](https://software.broadinstitute.org/software/igv/home)
* [ESEMBL Genome Browser](https://www.ensembl.org/index.html)

This tutorial was developed for use as part a series of workshops for neuroscience researchers, hence the example data and example genes are drawn from neuroscience field. However, the skills taught in this tutorial are applicable to all areas of research.

**Data:** [GTEX](https://gtexportal.org/home/) data as represented in the UCSC Genome Browser</br>
**Tools:** [UCSC genome Browser](https://genome.ucsc.edu/), [Integrative Genomics Viewer](https://software.broadinstitute.org/software/igv/home), [ESEMBL Genome Browser](https://www.ensembl.org/index.html)



-------------------------------
## Learning Objectives

At the end of this introductory workshop, you will :

* Understand the how different types of genomics and expression data are represented in Genome browsers
* know how to import and view your own data in a genome browser
* Know the main file types used for storing expression data BAM file type and
* Know how to index and sort BAM files using IGV tools
* Create, save and share custom views data in Genome browsers
*


-------------------------------
## Requirements and preparation

!!! attention "Important"</br>
    **Attendees are required to provide their own laptop computers.**  

If being delivered as a workshop, participants should install the software and data files below least one week before the workshop.  This should provide sufficient time for participants to liaise with their own IT support should they encounter any IT problems.  

Unless stated otherwise, recommended browsers are Firefox or Chrome (don't use Internet Explorer or Safari), and please ensure that your Browser is up to date.

Create a user account in UCSC genome browser
    Download and install [IGV](https://software.broadinstitute.org/software/igv/home) (free) on your laptop.

### Preparing your laptop prior to starting this workshop
1. Required software:
 * Ensure that ([Chrome](https://www.google.com/chrome/) or [FireFox](https://www.mozilla.org/en-US/) are installed and upto date)</br>
 * Download and install [IGV](https://software.broadinstitute.org/software/igv/download) (free) on your laptop.</br>
 * Check that the IGV software and data are correctly installed by executing this test(test)
3. Required data is downloaded as part of the tutorial exercises.


### Required Software
* [IGV](https://software.broadinstitute.org/software/igv/download)
* [Software x](https://...)

or

* No additional software needs to be installed for this workshop.

### Required Data
* [Data file](https://...)

or

* No additional data needs to be downloaded for this workshop.

-------------------------------
## Author Information
Written by: Victoria Perreau  
Melbourne Bioinformatics, University of Melbourne

Created/Reviewed: September 2020


-------------------------------
## Background

Little bit of history and context. Why is this important/useful...


-------------------------------

## Section 1: Intro to Genome Browsers

In this section you will ...

!!! attention "Important"
    **Please look at the [formatting template](/formatting_template/#formats-to-use) to see examples for formatting your documentation in line with Melbourne Bioinformatics training material e.g Questions and answers, code blocks, images etc.**  

Genome browsers are invaluable for viewing and interpreting many different types of data that can be anchored to genomic positions.  These include genomic variation, transcription, and the many types regulatory data such as methylation and transcription factor binding. In addition to being able to load your own data files the larger browsers also serve as data archives of valuable datasets facilitating visualisation and analysis of different data types of private and public data at the same time.  For example, coverage plots for cell type specific RNAseq data can be viewed at the same time as tissue level expression data and gene models illustrating transcript variants alongside genomic variation data. The UCCS genome browser Hub houses a number of important human brain gene expression datasets including GTEX and the LIBD developmental atlas.

Genome browsers allow viewing of one type of data in the context of other data which can reveal important information about gene regulation in normal development and disease. For example, they allow researchers to ask questions about whether a genomic variant could affect expression in specific transcript variants, cell types or tissues and contribute to hypothesis development relating to genotype phenotype relationships.  Note the vertical blue line in Figure below marking the position of a disease variant and showing that the two OMIM allelic variants known for BDNF are located in the coding exon, one of these variants also overlaps with the BDNF antisense transcript.

All researchers are therefore encouraged to become familiar with the use of at least one of the main browsers such as the UCSC ([UCSC Genome Browser] (https://genome.ucsc.edu/), RRID:SCR_005780), Ensembl  (Ensembl Genome Browser, RRID:SCR_013367), Epigenome browser at WashU (WashU Epigenome Browser, RRID:SCR_006208), IGV (Integrative Genomics Viewer, RRID:SCR_011793). They are designed for use by researchers without programming experience and have extensive tutorials and cases studies demonstrating the myriad of ways in which data can be loaded and interpreted to develop and support research hypothesis.

* Explore features of particular chromosomal regions
    * Investigate specific genes as well as collections of genes
    * Search for locations of sequences and markers
    * Retrieve annotation information for specific regions or genome-wide
    * View your own data in context of other annotations
    * Compare a region of one genome to genomes of other species

    ### Genome Build version number
    The Genome reference consortium
    https://www.ncbi.nlm.nih.gov/grc
    How does the nomenclature mean?
    https://genome.ucsc.edu/FAQ/FAQreleases.html

    For further info on Human Genome version updates I recommend you look at the updates and blog pages on the UCSC genome browser:
    https://genome.ucsc.edu/goldenPath/newsarch.html#2019
    http://genome.ucsc.edu/blog/patches/

-------------------------------
## Section 2: UCSC genome Browser

In this section we will become familar with e interface of the UCSC genome browser and explore a small number of tools and datasets available.

Weekly maintenance of the browser is at 5-6 pm pacific time = 11am-12pm Melbourne time. To ensure uninterrupted browser services for your research during UCSC server maintenance and power outages, bookmark a mirror site that replicates the UCSC genome browser.

## Section 3: ENSEMBL Genome Browser

In this section we will ...

## Section 4: IGV

In this section we will ...

Examine the new file by clicking on its <img src="../media/Galaxy-view.png" width=20 /> icon. We now have 2 columns instead of the 18 in the original file.

**1. Subheadding.**

* From the tool panel, click on **Text Manipulation -> Remove beginning** and set the following:
* "Remove First": *1*

Note the the new file is the same as the previous one without the header line.

**2. Make a histogram.**

* From the tool panel, click on **Graph/Display Data -> Histogram** and set the following:
* "Dataset": *Remove beginning on Data X*
* "Numerical column for X axis": *c2*


Click on the <img src="../media/Galaxy-view.png" width=20 /> icon of the histogram to have a look at it. Note there are a few peaks. Maybe these correspond to single, double and triple copy number of these contigs.

**3. Calculate summary statistics for contig coverage depth.**

* From the tool panel, click on **Statistics -> Summary Statisitics** and set the following:
* "Summary statistics on": *Remove beginning on Data X*



### Example 2: Convert Fastq to Fasta

This shows how to convert a fastq file to a fasta file. The tool creates a new file with the converted data.

**Converter tool**

* From the tool panel, click on **Convert Formats -> FASTQ to FASTA** and set the following:
* "FASTQ file to convert": *Typical Fastq File*
* Click **Execute**

This will have created a new Fasta file called FASTQ to FASTA on data 2.

### Example 3: Find Ribosomal RNA Features in a DNA Sequence


-------------------------------
## Additional reading
Links to additional recommended reading and suggestions for related tutorials.
