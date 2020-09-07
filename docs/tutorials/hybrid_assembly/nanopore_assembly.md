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

### Getting the data

1. **Make sure you have an instance of Galaxy ready to go.**
    * Navigate to the [Galaxy Australia server](https://usegalaxy.org.au/) and sign in if you have an account. 
2. **Copy an existing history**
    * The data you will need is available in an existing Galaxy history. You can create a copy of this history by clicking [here](https://usegalaxy.org.au/u/graceh1024/h/hybrid-de-novo-assembly) and using the import history '+' icon at the top right of the page. <br>
    <img src="../media/import_history.PNG" width="700">
3. **Look at the history you imported**
    * There are 4 files - Nanopore reads, a set of paired-end Illumina reads, and a reference genome for the organism we will assemble.
    * Will we use this reference genome to assess the quality of our assemblies and judge which methods work best. 

### Draft assembly with Flye + Nanopore reads
1. **Create the first draft assembly**<br><br>
We can use Flye to create an assembly from Nanopore reads.<br><br>
    * Making sure you are on the 'Analyse Data' tab of Galaxy, look for the tool search bar at the top of the left panel. 
    * Search for 'Flye' and select the tool 
    * We need to provide some information to Flye. Set the 'Input reads' parameter to nanopore_reads.fastq and 'estimated genome size' to 4m. Leave all else defualt. 
    * Run Flye by clicking 'execute' at the bottom of the page. 
    * Flye produces a number of outputs. We only need the 'consensus' fasta file. You can delete the other outputs. 
    * For clarity, the consensus draft assembly can be renamed to something which makes sense, like 'nanopore draft assembly'
<br><br>
2. **Assess draft assembly quality**<br><br>
We need to check if our assembly is good quality or not. It is paramount that genome assemblies are high-quality for them to be useful. <br> <br>
The supplied reference genome allows a direct comparison. We can use a tool call 'Quast' to compare our assembly to the reference genome.<br><br>
    * Search for the Quast tool in the tools panel. 
    * Parameters:
        - Contigs/scaffolds file = the nanopore draft assembly you just created
        - Use a reference genome? = Yes
        - Reference genome = reference_genome.fasta
        - All else default
    * Execute Quast by clicking 'execute' at the bottom of the page. 
    * We are mainly interested in one of the outputs - the HTML report
    * Open the report. It may look something like this:<br><br>
    <img src="../media/quast_draft_assembly.PNG" width="300">
    * Note the Genome fraction (%), # mismatches per 100 kbp, # indels per 100 kbp and # contigs information. 
<br><br>
We seem to have good coverage and not too many contigs, but our error rate is quite high.  We should be able to improve the error rate by polishing our assembly with high-accuracy Illumina reads by aligning them to our draft assembly and error correcting.
<br><br>
In this case we have a reference genome available, but this is not always the case. When our sample organism is unknown, we need another method to assess assembly quality. BUSCO analysis is one way to do this. 
<br><br>
BUSCO analysis uses the presence, absence, or fragmentation of key genes in an assembly to determine is quality. <br>
BUSCO genes are specifically selected for each taxonomic clade, and represent a group of genes which each organism in the clade is expected to possess. At higher clades, 'housekeeping genes' are the only members, while at more refined taxa such as order or family, lineage-specific genes can also be used. 
<br><br>
    * Find and select the Busco tool in the tools panel using the search bar.
    * We will assess our Nanopore draft assembly created by Flye.
    * In this tutorial, we will suspect that our organism is within the 'Bacillales' order. 
    * Parameters:
        * Sequences to analyse = our Nanopore draft assembly
        * Lineage = Bacillales
    * Leave all else default and execute the program.





        



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