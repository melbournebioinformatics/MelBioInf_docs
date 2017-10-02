# History creation instructions for Workflow tutorial

Use this set of instructions to create the base history for the "Extract Workflow" section of the workflows tutorial.

## Step 1: Import the raw datafiles

* Create a new blank history by clicking on the history menu <img src="../media/Galaxy-menu.png" width=20 />, then **Create New**
* Use the upload data tool to upload the data files from a remote repository..
    * Click **Get Data -> Upload File**
    * Click **Paste/Fetch data**
    * In the box paste the following two url's (one per line): *https://swift.rc.nectar.org.au:8888/v1/AUTH_377/public/COMP90014/Assignment1/bacterial_std_err_1.fastq.gz* and *https://swift.rc.nectar.org.au:8888/v1/AUTH_377/public/COMP90014/Assignment1/bacterial_std_err_2.fastq.gz*
    * Change the **Type** to *fastqsanger*
    * Click the **Paste/Fetch data** button again.
    * Paste *https://swift.rc.nectar.org.au:8888/v1/AUTH_377/public/COMP90014/Assignment1/Ecoli-O157_H7-Sakai-chr.fna*
    * Change the **Type** on this file to *fasta*
    * Click the **Start** button.
    * Click the **Close** button.
    
After the import is complete, you should have 3 files in your history. 2 fastq files and a fasta file.

## Step 2: Run BWA

Now we will run BWA on these files to map the reads to the reference.

* In the tools menu, click **NGS: Mapping -> Map with BWA**
* Set the following in the tool interface:
    * "Will you select a reference genome from your history or use a built-in index?": *Use a genome from history and build index*
    * "Use the following dataset as the reference sequence": *https://swift.rc.nectar.org.au:8888/v1/AUTH_377/public/COMP90014/Assignment1/Ecoli-O157_H7-Sakai-chr.fna*
    * "Select first set of reads": *https://swift.rc.nectar.org.au:8888/v1/AUTH_377/public/COMP90014/Assignment1/bacterial_std_err_1.fastq*
    * "Select second set of reads": *https://swift.rc.nectar.org.au:8888/v1/AUTH_377/public/COMP90014/Assignment1/bacterial_std_err_2.fastq*
* Click **Execute**

This will run BWA and result in an output compressed BAM file containing the mapping information for all of the reads versus the reference.

## Step 3: Run Freebayes

Now we will run Freebayes to call variants in out reads compared with the reference.

* In the tools menu, click **NGS: Variant Analysis -> Freebayes**
* Set the following in the tool interface:
    * "Choose the source for the reference genome": *History*
    * "BAM file": *Map with BWA on data 2, data 1, and data 3 (mapped reads in BAM format)*
    * "Use the following dataset as the reference sequence": *https://swift.rc.nectar.org.au:8888/v1/AUTH_377/public/COMP90014/Assignment1/Ecoli-O157_H7-Sakai-chr.fna*
* Click **Execute**

This will run Freebayes on your BAM file and will result in a variant calling format file (vcf). 

## Step 4: Filter the VCF file.

Now we will filter the vcf file to something sensible.

* In the tools menu, click **Filter and Sort -> Filter**
* Set the following in the tool interface:
    * "Filter": *FreeBayes on data 3 and data 4 (variants)*
    * "With following condition": *c6 > 500*
    * "Number of header lines to skip": *56*
* Click **Execute**

This will filter out VCF file.

You should now have the requisite history to enable you to complete the workflow extraction section of the workflows tutorial.

Return to it [here](galaxy-workflows.md)

    