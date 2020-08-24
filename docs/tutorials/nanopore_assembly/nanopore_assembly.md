<img src="../media/melbioinf_logo.png" width="350"> <img src="../media/PRIMARY_A_Vertical_Housed_RGB.png" width="150">

# Workshop title

Anticipated workshop duration when delivered to a group of participants is **2 hours**.  

For queries relating to this workshop, contact Melbourne Bioinformatics (bioinformatics-training@unimelb.edu.au).

## Overview

### Topic

* [x] Genomics
* [] Transcriptomics
* [ ] Proteomics
* [ ] Metabolomics
* [ ] Statistics and visualisation
* [ ] Structural Modelling
* [x] Basic skills


### Skill level

* [x] Beginner  
* [ ] Intermediate  
* [ ] Advanced  

This workshop is designed for participants with no command line knowledge. We will use a web-based platform called Galaxy to run all required  programs.


### Description

*Assemble a genome!<br />Learn how to create high-quality genome assemblies using the powerful combination of nanopore and illumina reads*

This tutorial uses "Flye" to perform de novo assembly of a Bacillus Subtilis genome using nanopore reads. "Pilon" will then be used to error correct our assembly and improve the accuracy of assembled contigs.  

**Data:** Nanopore reads, illlumina reads, Bacillus Subtilis reference genome<br />
**Pipeline:** Hybrid de novo genome assembly<br />
**Tools:** Flye, Pilon<br />


-------------------------------
## Learning Objectives

At the end of this introductory workshop, you will :

* Understand how nanopore and illumina reads can be used together to produce a high quality assembly
* Be familiar with genome assembly and polishing programs
* Be able to assemble an unknown, previously undocumented genome using nanopore and illumina reads! 


-------------------------------
## Requirements and preparation

**Attendees are required to bring their own laptop computers.**  

All data and tools are available on usegalaxy.org.au. You will need a computer to connect to and use their platform. 
Before the tutorial, navigate to https://usegalaxy.org.au/ and use your email to create an account. Click "Login or register" in the top navigation bar of galaxy to do this.  


### Preparing your laptop prior to starting this workshop
* No additional software needs to be installed for this workshop.


### Required Data
* No additional data needs to be downloaded for this workshop.


-------------------------------
## Author Information
Written by: Grace Hall  
Melbourne Bioinformatics, University of Melbourne

Created/Reviewed: August 2020


-------------------------------
## Background

Most modern bioinformatic analyses are not possible without high-quality reference genomes. They are the backbone that supports current and future discoveries, and our collection of complete genome sequences are one of humanities greatest achievements. The number of complete prokaryote genome sequences is growing rapidly due to their importance and relative ease of assembly. This information has enabled a giant step forward in our understanding of the world around us, and has created countless medical treatments to improve human health.

Most bacteria are unculturable. The genomes of these organisms must be assembled from DNA extracted from environmental or patient samples. We are usually unaware of which organisms are present in our samples, leaving de novo approaches our only viable option for assembly. By using nanopore reads as a scaffold, and illumina reads to error-correct and polish our assembly, we can reconstruct the genome of these organisms. Creating complete sequences of undocumented bacteria will reveal information about our interaction with the microbial world around and within us, and will lead to further medical discoveries so we may live happier, healthier lives. 


-------------------------------

## Section 1: Assembly with nanopore reads

In this section you will use Flye to create a draft genome assembly from nanopore reads. <br />
We will look at some features of our nanopore read set, perform assembly, then assess the quality of our assembly. 

Our process will be:

1. View nanopore read statistics with FastQC
2. Perform assembly using Flye
3. Look at assembly statistics
4. Check per-base accuracy by mapping contigs to known reference genome


-------------------------------
## Section 2: Title of section 2

In this section we will improve our assembly using illumina reads. 
Illumina reads can polish our assembly by:

* Correcting base errors in contigs
* Extending contigs
* Joining contig breaks
* Introducing new contigs


Our process will be:

1. View illumina read statistics with FastQC
2. Polish our assembly using Pilon
3. Look at the assembly statistics for our polished assembly
4. View the improvements to our assembly by mapping the polished contigs to the known reference genome

-------------------------------
## Additional reading
Links to additional recommended reading and suggestions for related tutorials.
None yet! 