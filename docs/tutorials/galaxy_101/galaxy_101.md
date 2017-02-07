<style>
  .image-header-gvl {
    -webkit-column-count: 3;
    -moz-column-count: 3;
    column-count: 3;

    -webkit-column-gap: 140px;
    -moz-column-gap: 140px;
    column-gap: 140px;
  }
</style>
<div class="image-header-gvl">
    <a href=http://genome.edu.au><img src="../media/gvl_logo.jpg" alt="GVL logo" width="112"/></a>
    <a href="http://galaxyproject.org"><img src="../media/GalaxyLogoHighRes.png" width="164" /></a>
    <a href="https://www.melbournebioinformatics.org.au/"><img src="../../../img/melbioinf_logo.png" alt="Melbourne Bioinformatics logo" width="164"/></a>
</div>



# Introduction to Galaxy


Written and maintained by [Simon Gladman](mailto:simon.gladman@unimelb.edu.au) - Melbourne Bioinformatics (formerly VLSCI)

## Background

Galaxy is a web based analysis and workflow platform designed for biologists to analyse their own data. It comes with most of the popular bioinformatics tools already installed and ready for use. There are many Galaxy servers around the world and some are tailored with specific toolsets and reference data for analysis of human genomics, microbial genomics, proteomics etc.

There are some introductory slides available [here](https://docs.google.com/presentation/d/1dzHagGkswjH7MOZ7OACVXGU-riBs33K3J5lWpnCpPhs/pub?start=false&loop=false&slide=id.g537781ff1_2_18).

Basically, the Galaxy interface is separated into 3 parts. The tool list on the left, the viewing pane in the middle and the analysis and data history on the right. We will be looking at all 3 parts in this tutorial.

<img src="../media/Galaxy_components.png" alt="Galaxy interface components diagram" width="70%" />

This workshop/tutorial will familiarize you with the Galaxy interface. It will cover the following topics:

* Logging in to the server
* Getting data into galaxy
* How to access the tools
* Using to use some common tools

-------------------------------

## Learning Objectives

At the end of this tutorial you should:

1. Be able to register on and login to a Galaxy server.
2. Be able to upload data to a Galaxy server from:
    * A file on your local computer
    * A file on a remote datastore with an accessible URL.
3. Be able use tools in Galaxy by:
    * Accessing the tool via the tool menu
    * Using the tool interface to run the particular tool
    * Viewing/accessing the tool output.

-------------------------------

## Section 1: Preparation.

The purpose of this section is to get you to log in to the server..

1. Go to the ip address of your GVL Galaxy server (or if you don’t have one, open [Galaxy-Tut](http://galaxy-tut.genome.edu.au) server) in Firefox or Chrome (your choice) - Please don’t use Internet Explorer or Safari.

2. If you have previously registered on this server just log in:

    * On the top menu select: **User -> Login**
    * Enter your password
    * Click **Submit**

3. If you haven’t registered on this server, you’ll need to now.

    * On the top menu select: **User -> Register**
    * Enter your email, choose a password, repeat it and add a (all lower case) one word name
    * Click **Submit**

---------------------------------------

## Section 2: Getting data into Galaxy

There are 2 main ways to get your data into Galaxy. We will use each of these methods for 3 files and then work use those 3 files for the rest of the workshop.

* Start a new history for this workshop. To do this:

    * Click on the history menu button (the <img src="../media/Galaxy-menu.png" width=20 /> icon) at the top of the Histories panel.
    * Select **Create New**

It is important to note that Galaxy has the concept of "File Type" built in. This means that each file stored needs to have its type described to Galaxy as it is being made available. Examples of file types are: text, fasta, fastq, vcf, GFF, Genbank, tabular etc.

We will tell Galaxy what type of file each one is as we upload it.

### Method 1: Upload a file from your own computer

With this method you can get most of the files on your own computer into Galaxy. (there is a size limit)

* Download the following file to your computer: *https://swift.rc.nectar.org.au:8888/v1/AUTH_377/public/galaxy101/Contig_stats.txt.gz*

    * From the Galaxy tool panel, click on **Get Data -> Upload File**
    * Click the **Choose File** button
    * Find and select the *Contig_stats.txt.gz* file you downloaded and click **Open**
    * Set the "file format" to *tabular*
    * Click the **Start** button
    * Once the progress bar reaches 100%, click the **Close** button

The file will now upload to your current history.

### Method 2: Upload a file from a URL

If a file exists on a web resource somewhere and you know its URL (Unique resource location - a web address) you can directly load it into Galaxy.

* From the tool panel, click on **Get Data -> Upload File**

    * Click on the **Paste/Fetch Data** button
    * Copy and paste the following web address into the URL/Text box: *https://swift.rc.nectar.org.au:8888/v1/AUTH_377/public/COMP90014/Assignment1/bacterial_std_err_1.fastq.gz*
    * Set the file format to *fastqsanger* (not fastqcsanger)
    * Click **Start**
    * Once the progress bar has reached 100%, click **Close**

Note that Galaxy is smart enough to recognize that this is a compressed file and so it will uncompress it as it loads it.

### Method 2 (again): Get data from a public URL

Now we are going to upload another file from the remote data source.

Repeat the above for: *https://swift.rc.nectar.org.au:8888/v1/AUTH_377/public/MRSA0252.fna*

Note that this file is a *fasta* file and not a *fastqsanger* file.

\showable{Reveal detailed instructions}{details}

* From the tool panel, click on **Get Data -> Upload File**

    * Click on the **Paste/Fetch Data** button
    * Copy and paste the following web address into the URL/Text box: https://swift.rc.nectar.org.au:8888/v1/AUTH_377/public/MRSA0252.fna
    * Set the file format to *fasta*.
    * Click **Start**
    * Once the progress bar has reached 100%, click **Close**

\endshowable

The DNA sequence of *Staphlococcus aureus MRSA252* will be loaded into your history as a fasta file.

Your history should now look like this.

<img src="../media/starting-history.png" width="300px" align="center" />

### The data

Though we aren't going to focus on the contents of these files and what they mean from a bioinformatics standpoint, here is a brief description of each one.

* *Contigs_stats.txt*

    * this file contains a table of summary data from a de novo genome assembly (the process of attempting to recover the full genome of an organism from the short read sequences produced by most DNA sequencing machines. )
    * The columns contain a lot of information but the ones we will be using indicate the amount of data (or coverage) that went into making up each piece of the final assembly.  

<!-- -->

* *bacterial_std_err_1.fastq.gz*

    * This file contains sequence reads as they would come off an Illunina sequencing machine. They are in [fastq](https://en.wikipedia.org/wiki/FASTQ_format) format.  

<!-- -->

* *MRSA0252.fna*

    * This file contains the genome sequence of *Staphylococcus aureus MRSA252*. It is in [fasta](https://en.wikipedia.org/wiki/FASTA_format) format.

---------------------------------------

## Section 3: Play with the tools

The purpose of this section is to get you used to using the available tools in Galaxy and point out some of the more basic manipulation tools.

Firstly however, you’ll notice that two of the files have very long and confusing names. So we might want to change them. To do this we need to “edit” the file. So:

1. Click on the <img src="../media/Galaxy-edit.png" width=20 /> icon (edit) next to the file in the history called: *https://swift.rc.nectar.org.au:8888/v1/AUTH_377/public/COMP90014/Assignment1/bacterial_std_err_1.fastq*
2. In the "Name" text box, give it a new name. Call it: *Typical Fastq File*
3. Click the **Save** button.

Repeat the process for the MRSA252 fasta file. Rename it to *MRSA252.fna*

Now that’s better. There was a lot of other functionality hidden behind that edit (<img src="../media/Galaxy-edit.png" width=20 />) icon. You can change a file’s data type, convert its format and many other things. Feel free to play around with them at a later date.

Ok, back to the tools..

### Example 1: Histogram and summary statistics

The first thing we are going to do is produce a histogram of contig read coverage depths and calculate the summary statistics from the Contig_stats.txt file. To do this we need to cut out a couple of columns, remove a line and then produce a histogram. This will introduce some of the text manipulation tools.

Click on the <img src="../media/Galaxy-view.png" width=20 /> icon of the *Contig_stats.txt* file to have a look at it. Note that there are 18 columns in this file. We want column 1 and column 6. To do this:

**1. Cut out column 1 and column 6.**

* From the tool panel, click on **Text Manipulation -> Cut** and set the following:
* Set "Cut Columns" to: *c1,c6*
* "Delimited by": *Tab*
* "Cut from": *Contig_stats.txt*
* Click **Execute**

Examine the new file by clicking on it’s <img src="../media/Galaxy-view.png" width=20 /> icon. We now have 2 columns instead of the 18 in the original file.

**2. Remove the Header lines of the new file.**

* From the tool panel, click on **Text Manipulation -> Remove beginning** and set the following:
* "Remove First": *1*
* "from": *Cut on data1*
* click **Execute**

Note the the new file is the same as the previous one without the header line.

**3. Make a histogram.**

* From the tool panel, click on **Graph/Display Data -> Histogram** and set the following:
* "Dataset": *Remove beginning on Data X*
* "Numerical column for X axis": *c2*
* "Number of breaks": *25*
* "Plot title": *Histogram of Contig Coverage*
* "Label for X axis": *Coverage depth*
* Click **Execute**

Click on the <img src="../media/Galaxy-view.png" width=20 /> icon of the histogram to have a look at it. Note there are a few peaks.. Maybe these correspond to single, double and triple copy number of these contigs.

**4. Calculate summary statistics for contig coverage depth.**

* From the tool panel, click on **Statistics -> Summary Statisitics** and set the following:
* "Summary statistics on": *Remove beginning on Data X*
* "Column or expression": *c2*
* Click **Execute**

You’ll note that the summary statistics tool failed and is red in the history. There was an error! If you click on the filename, and then the bug <img src="../media/Galaxy-bug.png" width=20 /> symbol, it will tell you what went wrong. (There is a missing python library.) At this point, you would normally contact your Galaxy server administrator.

### Example 2: Convert Fastq to Fasta

This shows how to convert a fastq file to a fasta file. The tool creates a new file with the converted data.

**Converter tool**

* From the tool panel, click on **Convert Formats -> FASTQ to FASTA** and set the following:
* "FASTQ file to convert": *Typical Fastq File*
* Click **Execute**

This will have created a new Fasta file called FASTQ to FASTA on data 2.

### Example 3: Find Ribosomal RNA Features in a DNA Sequence

This example shows how to use a tool called “barrnap” to search for rRNAs in a DNA sequence.

**1. Find all of the ribosomal RNA's in a sequence**

* From the tool panel, click on **Annotation -> barrnap** and set the following:
* "Fasta file": MRSA252.fna
* Click **Execute**

A new file called *barrnap on data 3* will be produced. It is a gff3 file. (This stands for genome feature format - version 3. It is a file format for describing features contained by a DNA sequence.) Change it’s name to something more appropriate (click on the <img src="../media/Galaxy-edit.png" width=20 /> icon.) There is also a STDERR output file from this tool - just ignore this one.

Now lets say you only want the lines of the file for the 23S rRNA annotations. We can do this using a Filter tool.

**2. Filter the annotations to get the 23S RNAs**

* From the tool panel, click on **Filter and Sort -> Select** and set the following:
* "Select lines from":  (whatever you called the barrnap gff3 output)
* "the pattern": *23S*     (this will look for all the lines in the file that contain “23S”)
* Click **Execute**

Now you have a gff3 file with just the 23S annotations!


## What now?

Remember how we started a new history at the beginning? If you want to see any of your old histories, click on the history menu button <img src="../media/Galaxy-menu.png" width=20 /> at the top of the histories panel and then select “Saved Histories.” This will give you a list of all the histories you have worked on in this Galaxy server.

That's it. You now know a bit about the Galaxy interface and how to load data, run tools and view their outputs. For more tutorials, see [http://genome.edu.au/learn](http://genome.edu.au/learn)
