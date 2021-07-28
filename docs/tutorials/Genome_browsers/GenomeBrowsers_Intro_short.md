![melbioinf_logo](../media/melbioinf_logo.png){: style="width:350px; padding-right:50px"}       ![unimelb_logo](../media/PRIMARY_A_Vertical_Housed_RGB.png){: style="width:150px"}

# Introduction to Genome Browsers

Anticipated workshop duration when delivered to a group of participants is **1 hour**.<br>
This is an abbreviated version of a [more extensive Genome Browser workshop](https://www.melbournebioinformatics.org.au/tutorials/tutorials/Genome_browsers/GenomeBrowsers_Intro/).  

For queries relating to this workshop, contact Melbourne Bioinformatics at:</br> bioinformatics-training@unimelb.edu.au.

## Overview
This tutorial will introduce you to the genome browser format and illustrate how some freely available genome browsers can be used to interrogate a variety of data types, such as gene expression, genomic variation, methylation and many more.

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


## Description

*Learn how to make the most of Genome Browsers !*

By focusing on gene expression, this hands on tutorial will provide beginners with an introduction to both the UCSC Genome browser and IGV (Integrated Genome Viewer).  Tools and public datasets will be used to illustrate how the expression of transcript variants can be investigated in different, tissues and cell types using public data, including human RNAseq data from GTEX and mouse cell type RNAseq data from Tabula Muris, as viewed within the UCSC genome browser. A subset of Single cell RNAseq data from the Allen Brain Atlas Celltax study will also be downloaded from SRA and visualised in IGV.
The data and genes used in this workshop are taken from the neuroscience field, however the analysis approaches and tools illustrated can be applied to many research areas.

***This tutorial is in three parts:***

* Section 1 Introduction to the general features of genome browsers. </br>
* Section 2 Hands on tutorial of the [UCSC Genome Browser](https://genome.ucsc.edu/)
* Section 3 Hands on tutorials of the [Integrative Genomics Viewer](https://software.broadinstitute.org/software/igv/home)


This tutorial was developed for use as part a series of workshops for neuroscience researchers, hence the data and example genes are drawn from neuroscience field and focused on analysis and visualisation  of expression data. However, the skills taught in this tutorial are applicable to all areas of research.

**Data:**

* [Single cell mouse cortex data from the Linnarsson lab](http://linnarssonlab.org/cortex/)
* [Tabular Muris](https://tabula-muris.ds.czbiohub.org/) data as represented in the UCSC Genome Browser
* [Celltax single cell expression atlas data](http://casestudies.brain-map.org/celltax) downloaded from [SRA](https://www.ncbi.nlm.nih.gov/sra)</br>

**Tools:**

* [UCSC genome Browser](https://genome.ucsc.edu/)
* [Integrative Genomics Viewer](https://software.broadinstitute.org/software/igv/home)

-------------------------------
## Learning Objectives

At the end of this introductory workshop, you will :

* Understand the how some types of genomics and expression data are represented in Genome browsers.
* Understand gene models and identify differences between transcripts variants.
* Examine the tissue/cell type expression profiles of a gene of interest in expression data.
* Know some basic files types used in Genome browsers and upload and view local BAM files.

-------------------------------
## Requirements and preparation

**Attendees are required to provide their own laptop computers.**  

If delivered as a workshop, participants should install the software and data files below prior to the workshop. Ensure that you provide sufficient time to liaise with your own IT support should you encounter any IT problems with installing software. Unless stated otherwise, recommended browsers are Firefox or Chrome.


### Preparing you and your laptop prior to starting this workshop
1. Recomended software: Download and install [IGV](https://software.broadinstitute.org/software/igv/download) (Free)
 * Ensure that ([Chrome](https://www.google.com/chrome/) or [FireFox](https://www.mozilla.org/en-US/) are installed and upto date)</br>
 * Create a user account in the [UCSC genome browser](https://genome.ucsc.edu/).
 * Required data is downloaded as part of the tutorial exercises.

### Required Data
* No additional data needs to be downloaded prior to this workshop.

-------------------------------
### Author(s) and review date
**Written by:** Victoria Perreau | Melbourne Bioinformatics, University of Melbourne.</br>

**Created:** October 2020</br>
**Reviewed and revised:** March 2021</br>

-------------------------------
## Genome Browser background

Genome browsers are invaluable for viewing and interpreting the many different types of data that can be anchored to genomic positions.  These include variation, transcription, the many types regulatory data such as methylation and transcription factor binding, and disease associations. The larger genome browsers serve as data archives for valuable public datasets facilitating visualisation and analysis of different data types. It is also possible to load your own data into some of the public genome browsers.

By enabling viewing of one type of data in the context of another, the use of Genome browsers can reveal important information about gene regulation in both normal development and disease, assist hypothesis development relating to genotype phenotype relationships.

All researchers are therefore encouraged to become familiar with the use of some of the main browsers such as:

* [The UCSC Genome Browser](https://genome.ucsc.edu/), (RRID:SCR_005780)
* [ESEMBL Genome Browser](https://www.ensembl.org/index.html), (RRID:SCR_013367)
* [Epigenome browser at WashU](https://epigenomegateway.wustl.edu/browser/), (RRID:SCR_006208)
* [Integrative Genomics Viewer (IGV)](http://software.broadinstitute.org/software/igv/), RRID:SCR_011793).

They are designed for use by researchers without programming experience and the developers often provide extensive tutorials and cases studies demonstrating the myriad of ways in which data can be loaded and interpreted to assist in develop and supporting your research hypothesis.

* [The UCSC Genome Browser Youtube channel](https://www.youtube.com/c/ucscgenomebrowser/videos)
* [Ensemble Browser webinar course](https://www.ebi.ac.uk/training/online/courses/ensembl-browser-series/)


Many large genomic projects also incorporate genome browsers into their web portals to enable users to easily search and view the data. These include:

* [GTEx](https://gtexportal.org/home/)
* [gnomAD](https://gnomad.broadinstitute.org/)

-------------------------------
## BDNF and TrkB signalling

This tutorial uses the a well known and important signalling pathway in the central nervous system (CNS) to illustrate some of the Genome browser tools and utility.

![TrkB-schema-eng](https://upload.wikimedia.org/wikipedia/commons/thumb/7/7a/TrkB-schema-eng.png/256px-TrkB-schema-eng.png){: align=left }  

Brain Derived Neurotrophic factor (BDNF) protein is an important neurotrophin responsible for regulating many aspects of growth and development in different cells within the CNS. TrkB is an important receptor that binds extracellular BDNF and propagates the intracellular signalling response via a tyrosine kinase. This TrkB receptor protein is encoded by the NTRK2 gene.  

The NRK2 gene expresses a number of different transcript variants in different cell types. The most well studied of these is the full length TrkB receptor referred to as TrkB, which is mainly expressed in neuronal cell types. The other transcript variants all express the same exons encoding the extracellular domain of the receptor (shown in the fugure here in green) but have truncated intracellular domains, which do not include the tyrosine kinase domain and thus activate different signalling pathways upon binding to BDNF.  None of these truncated protein products have been well studied, but the most highly expressed receptor variant is known as TrkB-T1, and is known to be highly expressed in astocytes.  

Since the transcript variants are differently expressed in different cell types within the CNS the NTRK2 gene is a very useful example for exploring cell type specific transcript expression in available public data.


**Major CNS cell types:**</br>

* Neuron  (yellow cell in the image below)
* Astrocyte  
* Oligodendrocyte  
* Microglia  
* Ependymal  


![1209 Glial Cells of the CNS-02](https://upload.wikimedia.org/wikipedia/commons/thumb/a/a4/1209_Glial_Cells_of_the_CNS-02.jpg/512px-1209_Glial_Cells_of_the_CNS-02.jpg)



-------------------------------

## Section 1: Introduction to Genome Browsers

Genome browsers rely on a common reference genome for each species in order to map data from different sources to the correct location. A consortium has agreed on a common numbering for each position on the genome for each species. However, this position will vary based on the version of the genome, as error correction and updates can change the numbering. Therefore it is very important to know which version of the genome your data of interest is aligned to.

The sequence for the human reference genome was accumulated up over many years from sequence data from many different sources and does not represent the sequence of one single person. Instead it is a composite of fragments of the genome from many different people.  Also, unlike the human genome which is diploid, the human genome is haploid.  That is there is only one copy of each chromosome. It therefore does not reflect the variation on the population, or even the most common variants in the human genome. Exploring variation within human genome is very important and facilitated by genome browsers but not covered in this workshop.

### Genome Build version number- further reading
* [The Genome reference consortium](https://www.ncbi.nlm.nih.gov/grc)</br>
* [What does the nomenclature mean?](https://genome.ucsc.edu/FAQ/FAQreleases.html)

For further info on Human Genome version updates I recommend you look at the updates and [blog pages on the UCSC genome browser](https://genome.ucsc.edu/goldenPath/newsarch.html#2019).

-------------------------------
## Section 2: The UCSC genome Browser interface

In this section we will become familiar with the web interface of the UCSC genome browser and explore some of the tools and public datasets available.

* Explore features of particular chromosomal regions
* Investigate specific genes as well as collections of genes
* Search for locations of sequences and markers
* Retrieve annotation information for specific regions or genome-wide
* View your own data in context of other annotations
* Compare a region of one genome to genomes of other species

Weekly maintenance of the browser is at 5-6 pm Thursdays Pacific time, which is equivalent to 11am-12pm AEST time. During this time the browser may be down for a few minutes. To ensure uninterrupted browser services for your research during UCSC server maintenance and power outages, bookmark one of the mirror sites that replicates the UCSC genome browser.

**Accessing the tools**:
Many of the tools that we will explore can be selected via multiple different routes within the browser interface. One way to access many tools is via from the top toolbar on a pull down list, other tools can be accessed from within the browser window.  In the following instructions a series of <ss>blue boxes</ss> is used to indicate successive lower levels from the pull down menu when starting with the top toolbar.  For example, the notation below indicates that you should select 'Genome Browser' from the top tool bar and then click on 'Reset all user settings'.

<ss>Toolbar</ss> <ss>Genome Browser</ss> <ss>Reset all user settings</ss>

**Accessing help and training**:
This workshop the UCSC genome browser is supported by rich training resource which has new material added regularly [youtube channel](https://www.youtube.com/channel/UCQnUJepyNOw0p8s2otX4RYQ/videos). To access training to further develop your skills and go to:
<ss>Toolbar</ss>  <ss>Help</ss> <ss>Training</ss>


### Getting started
1. **Open the Browser interface:**  
    * Navigate to the [UCSC genome Browser](https://genome.ucsc.edu/) and sign in if you have an account.<br>
<img src="../media/UCSC_Home.png" width="700">
2. **First reset the browser, so that we all see the same screen:**  
    <ss>Toolbar</ss> <ss>Genome Browser</ss> <ss>Reset all user settings</ss>
3. **Select and open the human Genome Hg38 at the default position, there are a few different ways to do this**  
    * <ss>Toolbar</ss> <ss>Genomes</ss> (this takes you to the Genome gateway page)  
        * Check that GRCh38 is selected in 'human assembly' and click on the blue <ss>GO</ss> box  
    * <ss>Toolbar</ss> <ss>Genomes</ss> <ss>Human GRCh38hg38</ss> (takes you directly to the genome)

    You should see this screen, opening at a position on the X chromosome of Human genome version GRCh38 showing the gene model for the ACE2 gene.
<img src="../media/UCSC_Hg38_opening.png" width="700">

4. **Familiarise yourself with the main areas of the interface and locate:**
    * The main Toolbar
    * Blue bar track collections (data of similar types are collected together under the same 'Blue bar' heading). Scroll down to see additional data collections and which ones are turned on as default.
    * Genome species and version number
    * Position box
    * Navigation tool buttons
    * Chromosome ideogram
    * Genome view window
    * Pre loaded tracks, track titles:
        * The grey bars on the left of the genome view can be used for selecting and configuring the tracks.
        * You can change the order of the tracks by dragging these grey bars up and down.
    * Turn tracks on and off:
        * You can hide tracks by right clicking on the grey bar or by turning them off in the Blue bar collections.
        * You have to click on a 'refresh' button to the changes to be reflected in the the genome view window.
    * View the **configuration page** specific to a track.  The configuration page gives you a lot of information about the data track and its colouring. You can open the configuration page for a track by:
        * clicking on the grey bar for the track or,
        * clicking on the track title in the Blue bar collection. More information and options is usually available by selecting the configuration page for a track via the track title in the Blue Bar collections.
    * Select white 'resize' button to fit the genome view window to your screen

5. **Customise your view by using the 'Configure' tool to change the font size to 12.** Use either method below to open the Configure tool.
    * <ss>Toolbar</ss> <ss>View</ss> <ss>Configure browser</ss> <ss>text size 12</ss> <ss>submit</ss>
    * or start by clicking on white <ss>configure</ss> button below the genome view window.

6. **Practice navigating around the genome view.**
    * move left and right both the navigation buttons and your mouse
    * zoom in and out using navigation buttons
    * zoom in to a region of interest using **'Drag-and-select'**:
        * using your mouse select a region of interest by clicking the ruler (position track) at the very top of the genome view window.
    * This is also how to access the **'highlight tool'** which you will use in a later exercise to highlight a region of interest.
        * Click on the down arrow next to the highlight colour to select a different colour.
<img src="../media/UCSC_drag-and-select.png" width="700">

###Understanding the gene models
<img src="../media/gene_model.png" width="700"><br>

#### NTRK2
First we are going to familiarise ourselves with the gene model representation of the different transcripts of NTRK2.  

1.  **Navigate to the NTRK2 gene position in GRCh38 and view the gene models**
    * You can navigate to a different region by typing in the position box.
        * If you know the specific location you are interested in type in the location using the format "chr#:1234-1234".
        * If you have a gene of interest you can type in the gene name (eg: NTRK2). Note the autocompleted suggestions that appear when  you start typing.  You can select from one of the suggestions or click <ss>go</ss> and select from a wider range of options.
        * Type (or copy and paste) **NTRK2** or **chr9:84,665,760-85,030,334** into the position box.
    * <ss>Hide all</ss> tracks by selecting the white button below the genome view.
    * Turn on only the 'Genecode v32 Genemodels' in 'full' viewing mode by selecting from the blue bar group labelled 'Genes and Gene predictions'.
    * Turn on 'Conservation' track to 'full'
        * Dont forget to click <ss>refresh</ss>.
    * When you have navigated to the NTRK2 gene, zoom out until you can view all of the 5' UTRs and 3' UTRs for all transcript variants for this gene. Then drag the view left and right to center (like in Google maps) or 'drag and select' the region to center the gene in the Genome view. You should see something like the image below.
<img src="../media/UCSC_NTRK2_GENCODE.png" width="700">

    * Which strand is the gene encoded on / transcribed from? (+ or - strand)
    * Identify the exons, introns and UTRs
        * *Do regions of conservation only occur were there are coding regions?*
        * *How many different transcripts variants are there for this gene?*
        * *How do they differ?*
    * Select a coding region (full height boxes) towards the 3'UTR of the gene.
        * zoom in to the region until you can see the letters of the amino acid sequence.
            * *Why are some amino acid boxes red or green?*
        * Zoom in again until you can see each amino acid number.  
            * *Why do different transcripts have different amino acid numbers?*
    * Note that one of the transcript names is in white text with a black background, this is the transcript you selected from the autocompleted list or the search results.
    * Change the **'view settings'** for the track. Switch between <ss>dense</ss> <ss>squish</ss> <ss>pack</ss> <ss>full</ss> to see how it changes the representation of the models.
        * Right click on the track grey bar in the left of the genome window to access view settings.
    * Go to the configuration  page for the Gencode v32 track and change the gene names to also reveal the 'Gencode transcript ID' in the label.
    * The transcript names are now too long to fit on the screen. Go to the genome view configuration page (like you did to chane the font size at the beginning of the workshop) and change the number of characters in the label so that you can see the entire transcript label.  

#### Gene model revision Quiz
Test your understanding of gene model representation by attempting this 6 questions in this [quiz](https://forms.gle/AnqR2igkm2xgEqEs6).

### Gene expression data
<img src="../media/coverage_plot.png" width="700">

1. **The FACS derived data from the [Tabular Muris](https://tabula-muris.ds.czbiohub.org/) cell type data is also availble in the UCSC genome browser and can be visualised as a coverage plot**
    * Start at the view of the NTRK2 gene in the human genome and navigate to the Ntrk2 gene in the mouse genome using the 'View in other genomes tool'.
        * <ss>Toolbar</ss> <ss>View</ss> <ss>In Other Genomes (Convert)</ss>
        * select New Genome:Mouse , New Assembly:GRC38/mm10, click on 'Submit'
        * select the region with the greatest homology
    * Configure the Tabular Muris track by selecting it from the blue bar collection.
        * Hide 'Cell expression'
        * Select 'Genome coverage' to full
        * Select 'submit'<br>
**This can look like a bit too much data to manage as there are very many tracks and the default track height is set very high. But its easy to simplify it by filtering to show only a few cell types of interest.**<br>
    * Right click on the grey bar to 'configure the track set'.
        * Change track height to 30
        * for 'data view scaling' select group autoscale
        * clear all the subtracks and then manually select only a few cell types of interest:
            * astrocyte Cv
            * Bergmann glial Cv
            * microglia Cv
            * neuron Cv
            * oligodendrocyte Cv
            * OPC Cv
    * Which cell type has the highest level expression in this dataset?
    * Change the 'Data view scaling' to autoscale to dataview.
        * Export a PDF image of the genome view: <ss>Toolbar</ss> <ss>View</ss> <ss>PDF/PS</ss> select 'Download the current browser graphic in PDF'
    <img src="../media/mm10_Ntrk2.png" width="600">
    * Which cell type(s) express the long and short transcripts for NTRK2?

2. **We are going to examine the expression of transcript variants of NTRK2 in different cell types in the mouse brain cortex [Linnarsson lab](http://linnarssonlab.org/cortex/)**
    * The data that is publicly available for viewing in the UCCS genome browser but is not housed in the UCSC genome browser. You must first access it from the the Linnarsson lab data page.
    * This RNAseq data is stranded, meaning you can see if the transcript data is from the + or - strand.
    * Go to the **[Public data page](http://linnarssonlab.org/cortex/)** where you can search for cell expression profiles for individual genes.
    * Click on the 'Browse the genome' blue text near the bottom of the page.
    * This loads 18 different tracks, one for each cell type. The default setting for expression range is quite high and most gene expression cannot be visualised with these settings. Each track must be also be configured individually rather than as a group, which I will demonstrate but it takes a lot of time.
    * I have created a version of this data where each track autoscales which can make it quicker to determine what expression range would be ideal for visualising the expression of an individual gene.  You can access this custom track set at this **[link](https://genome.ucsc.edu/cgi-bin/hgTracks?db=mm10&lastVirtModeType=exonMostly&lastVirtModeExtraState=&emGeneTable=knownGene&virtModeType=exonMostly&virtMode=1&nonVirtPosition=chr13%3A58807697%2D59133970&position=chr13%3A58807697%2D59133970&hgsid=934998079_avsgVoOzSJintErsCjEIcwdiQOQc)**.
    * For an exercise, select 2 or three cell types and adjust the scale to best reflect differences in gene expression of Ntrk2 between these cells. Save this session and share it.



## Section 3: IGV

In this section we will download a BAM file of gene expression data from SRA and view it in the [Integrated Genome Viewer (IGV)](https://igv.org/).  

For simplicity, and in order to save time if you have not yet installed IGV on your computer, you can view the BAM files  the IVG web browser application. However you cannot generate the Sashimi plot view in the web browser version of IGV.

BAM files must first be sorted and indexed before they can be loaded into genome viewers and IGV has tools to do this in the desktop app without having to use command line. However these tools are not available in the [IGV web app](https://igv.org/app/). Therefore I have also provided the sorted BAM and index files for download, but I encourage you to download them yourself from SRA and sort and index them using the instructions below when you have time.<br>
**[BAM and Bai file download link](https://zenodo.org/record/4636667#.YFyAuGQzZqt)**

The expression data we are using for this exercise is from the mouse [Celltax single cell expression atlas](http://casestudies.brain-map.org/celltax) published by the Allen Brain Institute.  The cell tax vignette has an expression browser that displays gene level expression as a heat map for any gene of interest,  The readsets (fastq files) and aligned data (BAM files) for 1809 runs on single cells are also available for down load from SRA.

The SRA study ID for this study is [SRP061902](https://www.ncbi.nlm.nih.gov/sra?term=SRP061902) and individual runs from this study are easily selected by viewing the samples in the 'RunSelector', if you wish to identify particular cell types of interest.  For this exercise, I have already identified a few samples that we will download in order to illustrate navigating data in IGV by looking at the expression of NTRK2 in the same cell types we have discussed in earlier exercises.<br>

For each cell type we will down load a .BAM file containing only the reads from a single chromosome of interest. Using a reduced dataset for demonstrations cuts down on data transfer and processing time.

For each SRA run in the table below open the link in a new tab to down load the data. Not many readsets in SRA have aligned data available for down load but this data set does.

| Cell type  | SRA run  | Vignette Cell ID|
| :--------- | :------- |:------- |
| **astrocyte**   | [SRR2138962](https://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?run=SRR2138962)  | D1319_V |
| **astrocyte**   | [SRR2139935](https://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?run=SRR2139935)  | A1643_VL |
| **neuron** | [SRR2139989](https://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?run=SRR2139989)  | S467_V4 |
| **neuron** | [SRR2140047](https://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?run=SRR2140047) | S1282_V |

1. **Download BAM files from SRA**
    * Click on the 'Alignment' tab
    * Note that the data is aligned to the mouse GRCm38 genome (mm10).
    * Select the chromosome of interest. For NTRK2 in mouse it is chr13
    * For 'Output this run in:' select **BAM** and click on 'format to:' **File**
    * Rename the downloaded file to include the cell type, to avoid confusion. eg: **SRR2138661_astrocyte_chr13.bam**

2. **Use IGV tools to SORT and INDEX the BAM files** Store sorted BAM files and index files in the same folder.
    * Open IGV and select Tools / Run igvtools... from the pull town menus.
    * Select **'Sort'** from the Command options and use the brows options to select the BAM file you just downloaded and click **'Run'**
    * Without closing the igvtools window now select the command **'Index'** and Browse to find the BAM file you just sorted.  It will have the same file name with 'sorted' added to the end. eg SRR2138661_astrocyte_chr13.sorted.bam
    * The resulting index file will have the file name : SRR2138661_astrocyte_chr13.sorted.bam.bai
        * It is essential that the index file for a BAM file has the same name and is located in the same folder as its BAM file. If not the IVG software will not be able to open the BAM file.<br>

<img src="../media/BAM_bai_files.png" width="300">

3. **View the BAM files in IGV**
    * Select the Mouse (mm10) genome from the genome box in the top right hand corner.
    * Select File / Load from File... and select all 4 '_chr13.sorted.bam' files only (use command to select more then one file at a time).
    * select open - but don't expect to see any data yet. The genome view window opens on a whole chromosome view as default but it wont show any data until the view region is small enough to show all data in the current view.
    * Type the gene name 'NTRK2' into the search window.
    * Expand the Refseq gene model track by right clicking it to see all the splice variants
    * The gene and thus the genome view is 328kb and the default setting for viewing data is only 100kb.  So unless you have already changed your settings alignment data will not get be showing.
    * zoom into the region of a coding exon by selecting in the numbered location track at the top of the genome view.
    * To see the whole gene in the genome window at the same time you may need to change the preferences.  
    * Go to View / Preferences and select the 'Alignments tab'.
    * Change the visibility range threshold to 400kb.

<img src="../media/IGV_preferences.png" width="400">

**You may need to change this back to a smaller range in the future if you are working with large datasets and/or small amounts of memory on your computer.**

<img src="../media/igv_Ntrk2.png" width="700">

4. **Export images**
    * The Genome view above can be exported by selecting 'File / Save image...'  from the tool bar.
    * To export the Sashimi plot below:
        * Right click on one of the junction tracks and select 'Sashimi Plot' from the poll down menu.
        * Select the tracks you want in your final image.
        * There are some data filtering and style adjustments you can make to the Sashimi plot. Right click on each track to access the menu options. Some changes apply to each track individually and some to all tracks.
            *
<img src="../media/Sashimi.png" width="700">


-------------------------------
## Additional reading

**IGV**
https://rockefelleruniversity.github.io/IGV_course/presentations/singlepage/IGV.html
