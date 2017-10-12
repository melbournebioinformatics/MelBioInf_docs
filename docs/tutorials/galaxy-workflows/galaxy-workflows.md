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

# Galaxy Workflows


Written and maintained by [Simon Gladman](mailto:simon.gladman@unimelb.edu.au) - Melbourne Bioinformatics (formerly VLSCI)

## Background

This workshop/tutorial will familiarise you with the Galaxy workflow engine. It will cover the following topics:

* Logging in to the server
* How to construct and use a workflow by various methods
* How to share a workflow

---------

## Section 1: Preparation.

The purpose of this section is to get you to log in to the server.. For this workshop, you can use a Galaxy server that you have created by following the steps outlined in: [Launch a GVL Galaxy Instance]() or you can use the public [Galaxy-mel](http://galaxy-mel.genome.edu.au)


**Go to Galaxy URL of your server in Firefox or Chrome (your choice, please don't use IE or Safari.)**

* If you have previously registered on this server just log in:
    * On the top menu select: **User -> Login**
    * Enter your password
    * Click the **Submit** button

<!-- -->

* If you haven’t registered on this server, you’ll need to now.
    * On the top menu select: **User -> Register**
    * Enter your email, choose a password, repeat it and add a (all lower case) one word name
    * Click on the **Submit** button.

-----------

## Section 2: Create and run a workflow.

This section will show you two different methods to create a workflow and then how to run one.

### Import the workflow history

In this step we will import a shared history to our workspace so we can extract a workflow from it. This will only work on Galaxy servers which have the history available on it. If yours doesn't have the appropriate history, there are instructions to create it [here](history_creation.md). Alternatively, you can extract a workflow from any history you have in your "Saved Histories."

* From the menu at the top of the Galaxy window, click **Shared Data -> Histories**
* Find the history called "*workflow_finished*" and click on it.
* Then click on **Import history** at the top right of the screen.
* Change the name if you wish and then click **Import**

This history should now be in your history pane on the right.


<!--


This step will show you how to import a history from a remote source into your own workspace. We will be using this history to build a workflow.

* From the histories menu <img src="../media/Galaxy-menu.png" width=20 />, click **Import from File**.
* Enter *https://swift.rc.nectar.org.au:8888/v1/AUTH_377/public/Workflow-finished-history.tar.gz* in the URL box.
* Click **Submit**

You might have to wait for a bit, then:

* From the history menu <img src="../media/Galaxy-menu.png" width=20 />, click **Saved Histories**
* Select the history: *Imported: Workflow-finished*


-->
## Workflow creation: Method 1

We will create a workflow from an existing history. You can use this method to make a re-useable analysis from one you’ve already done. i.e. You can perform the analysis once and then create a workflow out of it to re-use it on more/new data. We will create a workflow from the history you imported in step 1 above. The footnote below explains the steps that created this history. These are the steps we will mimic in the workflow.

Make sure your current history is the one you imported in Section 2 - Step 1 above (*imported: workflow_finished.*) If not, switch to it. If you couldn't import it, you should be able to complete this step with any suitable history.

**Now we will create the workflow.**

* Click on the histories menu button <img src="../media/Galaxy-menu.png" width=20 /> at the top of the history pane.
* Click **Extract Workflow**

You will now be shown a page which contains the steps used to create the history you are extracting from. We use this page to say what to include in the workflow. We want everything here so we’ll just accept the defaults and:

* Change the Workflow name to something sensible like “Basic Variant Calling Workflow”
* Click **Create Workflow**

The workflow is now accessible via the bottom of the tool pane by clicking on All Workflows.


### Some discussion

Have a look at your workflow. Click on its button in the workflow list. It’s tool interface will appear. You can now run this workflow any time you like with different input datasets. NOTE: The input data sets must be of the same types as the original ones. i.e. In our case, two fastq reads files and one fasta reference sequence.

More interesting though is to:

* Click on the  **Workflows** link in the top menu
* Click on the down arrow on your workflow’s button.
* Click **Edit**

A visualisation of your workflow will appear. Note the connections and the steps.

Next we’ll go through how to create this workflow using the editor..

## Workflow Creation: Method 2

We will now create the same read mapping/variant calling workflow using the editor directly. The workflow needs to take in some reads and a reference, map the reads to the reference using BWA, run Freebayes on the BAM output from BWA to call the variants, and finally filter the resulting vcf file.

### Step 1: Create a workflow name and edit space.

* Click on **Workflow** in Galaxy’s menu.
* Click on the **Create New Workflow** button.
* In the "Workflow Name" text box type: *Variants from scratch*
* Click the **Create** button.

You should now be presented with a blank workflow editing grid.

### Step 2: Open the editor and place component tools

**Add three input datafiles.**

* In the Workflow control section of the tool pane, click on **Inputs -> Input dataset** three times.
* Spread them out towards the left hand side of the workflow grid by clicking and dragging them around.
* For each one, change their name.
    * Click on each input box in turn
    * In the right hand pane (where the history usually is), change the name to:
        * Reference data
        * Reads 1
        * Reads 2 - respectively.
    * If you click on **Input dataset** in the tan box at the top of the right hand panel on the screen, you can change the name of the box as it appears on the editing layout. You need to give each one a more sensible name and press *Enter* to make the change.

**Add in the BWA mapping step.**

* Click on **NGS: Mapping -> Map with BWA** in the tool pane. BWA will be added to the workflow grid.
* In the right hand pane (where the history usually is), change the following parameters.
    * Change “Will you select a reference genome from your history or use a built-in index?:” to *Use a genome from history.*
    * Note that the BWA box on the grid changes to match these settings.

**Connect the tools together.**

* Click and drag on the output of one of the input dataset tools to each of the input spots in the BWA tool. (Make connections.)
    * "Reference data" output to reference to "Use the following dataset as the reference sequence"
    * "Reads 1" output to "Select first set of reads"
    * "Reads 2" output to "Select second set of reads"

**Add in the Freebayes (Variant Calling step.)**

* Click on **NGS: Variant Calling -> Freebayes**
* In the right hand pane, change the following:
    * "Choose the source for the reference list:"  to *History*
    * Connect "Map with BWA" bam output to "Freebayes’" bam input.
    * Connect the "Reference data" output to "Freebayes’" *Use the following dataset as the reference sequence* input.

If you’re keen - **Note: this is optional.** Also change the following parameters the right hand pane for freebayes to make it a bit more sensible for variant calling in bacterial genomes.

* "Choose parameter selection level": *Complete list of all options*
* "Population model options": *Set population model options*
* "Set ploidy for the analysis": *1*
* "Input filters": *Set input filters*
* "Exclude alignments from analysis if they have a mapping quality less than": *20*
* "Exclude alleles from analysis if their supporting base quality is less than": *20*
* "Require at least this fraction of observations … to evaluate the position": *0.9*
* "Require at least this count of observations .. to evaluate the position": *10*
* "Population and mappability priors": *Set population and mappability priors*
* "Disable incorporation of prior expectations about observations": *Yes*

**Add in the Filter step.**

* Click on **Filter and Sort - > Filter**
* Connect "Freebayes’" output_vcf to the "filter" input.
* In the right hand pane, change the following:
    * "With the following condition": *c6 > 500*
    * "Number of header lines to skip": *56*


Phew! We’re nearly done! The only thing left is to select which workflow outputs we want to keep in our history. Next to each output for every tool is a star. Clicking on the stars will select those files as workflow outputs, everything else will be hidden in the history. In this case we only really want the BAM file and the final variants (vcf file.) Therefore:

**Select workflow outputs.**

* Click on the star next to "Map with BWA’s" bam file output.
* Click on the star next to "Filter’s" output vcf.

### Step 3: Save it!

Click on the <img src="../media/Galaxy-menu.png" width=20 /> at the top of the workflow grid and select **Save**.

Congratulations. You’ve just created a Galaxy workflow.

Now to run it!

## Running the workflow

We will now make a new history called "Test" and run the workflow on it’s data.

### Create the new history

* From the Histories Menu, select **Copy Datasets**
* Select the 2 x fastq files and the Ecoli .fna file.
* Under destination history, enter *Test* into "New History Named:"
* Click **Copy History Items** button
* Click on the link to the new history in the green bar at the top of the screen

### Run the workflow

* On the tools pane, click **All Workflows**
* Select the *Variants from scratch* workflow (or whatever you called it.)
* Give it the correct files.
* *Ecoli ... .fna* for Reference data
* *bacterial_std_err_1.fastq* for Reads 1
* *bacterial_std_err_2.fastq* for Reads 2
* Click **Run Workflow**

Your workflow will now run. It will send the right files to the right tools at the right time to the cluster (compute engine on your machine) and wait for them to finish. Watch as they turn yellow then green in turn.

## What now?

Where to start? There’s so much you can do with Workflows. You can even run them on multiple file sets from the one setup.
