<img src="../media/melbioinf_logo.png" width="350"> <img src="../media/PRIMARY_A_Vertical_Housed_RGB.png" width="150">

# Hybrid genome assembly - nanopore and illumina

Anticipated workshop duration when delivered to a group of participants is **2 hours**.  

For queries relating to this workshop, contact Melbourne Bioinformatics (bioinformatics-training@unimelb.edu.au).

## Overview

### Topic

* [x] Genomics
* [ ] Transcriptomics
* [ ] Proteomics
* [ ] Metabolomics
* [ ] Statistics and visualisation
* [ ] Structural Modelling
* [x] Basic skills


### Skill level

* [x] Beginner  
* [ ] Intermediate  
* [ ] Advanced  

This workshop is designed for participants with no command line knowledge. A web-based platform called Galaxy will be used to run our analysis.


### Description

*Assemble a genome!<br>Learn how to create and assess the quality of high-quality genome assemblies using the powerful combination of nanopore and illumina reads*

This tutorial explores how long and short read data can be combined to assemble bacterial genomes, and how the quality of these assemblies can be assessed. Termed 'hybrid assembly', this method is now the gold-standard for creating high-quality genome assemblies from scratch. We will perform de-novo (assembly from scratch, no reference genome) from Nanopore and Illumina reads of a single bacterial organism. Two hybrid assembly methods will be performed, and the quality of resulting assemblies will be assessed.

**Data:** Nanopore reads, Illlumina reads, bacterial organism (Bacillus subtilis) reference genome<br>
**Tools:** Flye, Pilon, Unicycler, Quast, BUSCO<br>
**Pipeline:** Hybrid de novo genome assembly - Nanopore draft Illumina polishing<br>
**Pipeline:** Hybrid de novo genome assembly - Unicycler<br>



-------------------------------
## Learning Objectives

At the end of this introductory workshop, you will :

* Understand how Nanopore and Illumina reads can be used together to produce a high quality assembly
* Be familiar with genome assembly and polishing programs
* Learn how to assess the quality of a genome assembly, regardless of whether a reference genome is present or absent
* Be able to assemble an unknown, previously undocumented genome to high-quality using Nanopore and Illumina reads! 


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

Created/Reviewed: September 2020


-------------------------------
## Background

Most modern bioinformatic analyses are not possible without high-quality reference genomes. They are the backbone that supports current and future discoveries, and our collection of complete genome sequences are one of humanity's greatest achievements. The number of complete prokaryote genome sequences is growing rapidly due to their importance and relative ease of assembly. This information has enabled a giant step forward in our understanding of the world around us, and has created countless medical treatments to improve human health.

Most bacteria are unculturable. The genomes of these organisms must be assembled from DNA extracted from environmental or patient samples. We are usually unaware of which organisms are present in our samples, leaving de novo approaches our only viable option for assembly. Using the powerful combination of Nanopore and Illumina reads, we can reconstruct the genome of these organisms. Creating complete genome assemblies of undocumented bacteria will reveal information about our interaction with the microbial world around and within us, and will lead to further medical discoveries so that we may live happier, healthier lives. 


-------------------------------

## Section 1: Nanopore draft assembly, Illumina polishing

In this section you will use Flye to create a draft genome assembly from nanopore reads. <br />
We will look at some features of our nanopore read set, perform assembly, then assess the quality of our assembly. 

Our process will be:

1. Perform draft assembly using Nanopore reads + Flye
2. Assess assembly quality using Quast (compare to reference genome) and BUSCO analysis
3. Polish and improve our assembly with Illumina reads + Pilon
4. Reassess assembly quality and compare


-------------------------------
## Section 2: Purpose-built hybrid assembly tool - Unicycler

In this section we will use a purpose-built tool called Unicycler to perform hybrid assembly.

Our process will be:

1. Perform assembly with Nanopore reads, Illumina reads + Unicycler
2. Assess assembly quality using Quast (compare to reference genome) and BUSCO analysis
3. Compare our two approaches


-------------------------------
## Additional reading
Links to additional recommended reading and suggestions for related tutorials.<br>
Flye: https://github.com/fenderglass/Flye/blob/flye/docs/USAGE.md#algorithm<br>
Pilon: https://github.com/broadinstitute/pilon/wiki/Methods-of-Operation<br>
Unicycler: https://github.com/rrwick/Unicycler<br>
Quast: https://academic.oup.com/bioinformatics/article/29/8/1072/228832<br>
BUSCO analysis: https://academic.oup.com/bioinformatics/article/31/19/3210/211866<br>