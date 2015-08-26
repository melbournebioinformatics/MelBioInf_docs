# RNA-Seq Experimental Design

## What is RNA-seq?

RNA-seq is a method of measuring gene expression using shotgun sequencing. The
process involves reverse transcribing RNA into cDNA, then sequencing fragments
on a high-throughput platform such as Illumina to obtain a large number of
short reads. For each sample, the reads are then aligned to a genome, and the
number of reads aligned to each gene or feature is recorded.

A typical RNA-seq experiment aims to find differentially expressed genes
between two conditions (e.g. up and down-regulated genes in knock-out mice
compared to wild-type mice). RNA-seq can also be used to discover new
transcripts, splice variants, and fusion genes.


## Why is a good experimental design vital?

An RNA-seq experiment produces high dimensional data. This means we get a huge
number of observations for a small number of samples. For example, the
expression of ~20,000 genes could be measured for 6 samples (3 knock-out and 3
wild-type). A frequently used approach to analyse RNA-seq data is to fit each
gene to a linear model where for each of the 20,000 genes, parameters need to
be estimated using a small number of observations. To complicate matters, each
measurement of gene expression is comprised of a mix of biological signal and
unwanted noise. Thus, in order to perform a robust statistical analysis, the
methodology must be carefully designed.

Before you begin any RNA-seq experiment, some questions you should ask
yourself are:

 - Why do you expect to find differentially expressed genes in the particular
   tissue?
 - What types of genes do you expect to find differentially expressed?
 - What are the sources of variability from your samples?
 - Where do you expect most of your variation to come from?

A coherent experimental design is the groundwork of a successful experiment.
You should invest time and thought in designing a robust experiment as failing
to think this step through can lead to unusable data and wasted time, money,
and effort.

It is also useful to think about the statistical methods you will use to
analyse the data. If you're planning to bring a data analyst or
bioinformatician onboard for data analysis, you should include him or her in
the experimental design stage.


## Terminology

Before progressing, it may be useful to define some terms which are commonly
used in RNA-seq.

<table>
  <tr>
    <td>Variability:
    <td>A measure of how much the data is spread around. Variance is
    mathematically defined as the average of the squared difference between
    observations and the expected value. Simply put, a larger variance means
    it is harder to identify differentially expressed genes.
  <tr>
    <td>Feature:
    <td>A defined genomic region. Usually a gene in RNA-seq, but can also
    refer to any region such as an exon or an isoform. In RNA-seq, an estimate
    of abundance is obtained for each feature.
  <tr>
    <td>Biological replicates:
    <td>Samples that have been obtained from biologically separate samples.
    This can mean different individual organisms (e.g. tissue samples from
    different mice), different samplings of the same tumour, or different
    population of cells grown separately from each other but originating from
    the same cell-line. For example, the samples obtained from three different
    knock-out mice could be considered biological replicates in a knock-out
    versus wild-type experiment. A biological replicate combines both technical
    and biological variability as it is also an independent case of all the
    technical steps.
  <tr>
    <td>Technical replicates:
    <td>Samples in which the starting biological sample is the same, but the
    replicates are processed separately. For example, if a biological sample
    is divided and two different library preps are processed and sequenced,
    those two samples would be considered technical replicates.
  <tr>
    <td>Covariate:
    <td>The term 'covariate' is often used interchangeably with 'factor' or
    'variable' in RNA-seq. The term refers to a property of the sample which
    may have some influence on gene expression and should be represented in
    the RNA-seq model. Covariates in RNA-seq are often categorical (e.g.
    treatment condition, sex, batch), but continuous factors are also possible
    (e.g. time points, age). A linear model will contain terms to represent
    the relationships between covariates and each sample. Each possible value
    a factor can take is called a level (e.g. 'male' and 'female' are two
    levels in the factor 'sex'). Factors can either be directly of interest
    to the experiment (e.g. treatment condition) or not of interest (also
    known as nuisance variables) (e.g. sex, batch). The purpose of covariates
    is to explain the variance seen in samples.
  <tr>
    <td>Confounding variable:
    <td>A confounding variable is a nuisance variable that is associated with
    the factor of interest. Possible confounding factors should be controlled
    for so they don't interfere with analysis. For example, if all knock-out
    mice samples were harvested in the morning and all wild-type mice samples
    were harvested in the afternoon, the time of sample collection would be a
    confounding factor as the effects from sample collection time and from the
    knock-out cannot be separated.
  <tr>
    <td>Statistical power:
    <td>The ability to identify differentially expressed genes when there
    really is a difference. This is partly dependent on variance and therefore
    is affected by the number of replicates available and sequencing depth.
</table>



## The importance of replicates to estimate variance

When performing a differential gene expression analysis, we look at the
expression values of each gene and try to determine if the expression is
significantly different in the different conditions (e.g. knock-out and
wild-type). The ability to distinguish whether a gene is differentially
expressed is partly determined by the estimates of variability obtained by
using multiple observations in each condition.

Variability is present in two forms: technical variability and biological
variability.

Combined biological and technical variability is measured using biological 
replicates. Biological variability is the main source of variability and is 
due to natural variation in the population and within cells. This includes 
different individuals having different levels of a particular gene and the 
stochastic nature of expression levels in different cells.

Technical variability is measured using technical replicates. Technical
variability is often very small compared to biological variability.  Usually
the question is whether an observed difference is greater than the total 
variability (i.e. significant).  As combined variability is measured by 
biological replicates technical replicates are only important if you need to
know the degree of biological variability or technical variability.  An
example of wanting technical variability would be method development. The
main source of technical variation comes from RNA processing and from 
library prep. Variability from sequencing in different flow cells or different
lanes is usually minimal. Generally, creating technical replicates from multiple
library preps is unnecessary for RNA-seq experiments.

The amount of variance between your biological replicates will affect the
outcome of your analysis. Ideally, you aim to have minimal variability between
samples so you only measure the effect of the condition of interest. Too much
variability between samples can drown out the signal of truly differentially
expressed genes. Controlling for possible confounding factors between
conditions is also important to prevent falsely attributing differential
expression to the condition of interest.

Strategies to minimise variation between samples and to control confounding
variables include:

 - choosing organisms from the same litter,
 - choosing organisms of the same sex if possible,
 - using a constant sample collection time,
 - having the same laboratory technician perform each library prep,
 - randomising samples to prevent a confounding batch effect if all samples
   can't be processed at one time.

If variation between samples can not be removed it should be balanced between
conditions of interest as much as possible, and carefully recorded to allow
its effect to be measured and potentially removed during analysis.


## How many replicates and how many reads do I need?

Two very common question asked are:

 - how many biological replicates do I need, and
 - what sequencing depth is needed for each sample

in order to have enough statistical power for my RNA-seq experiment?

These questions cannot be precisely answered without a pilot study. A small
amount of data (minimum of two biological replicates for each condition with
at least 10M reads) can estimate the amount of biological variation, which
determines how many biological replicates are required. Performing a pilot
study is highly recommended to estimate statistical power and identify
possible problems before investing more time and money into the project.

[Scotty](http://euler.bc.edu/marthlab/scotty/scotty.php) is a web-based tool
that uses data generated from a pilot study to optimize a design for
statistical power. With a limited budget, one must balance sequence coverage
and number of biological replicates. Scotty also has a cost estimate feature
which returns the most powerful design within budget constraints.

As a general rule, the number of biological replicates should never be below 3.
For a basic RNA-seq differential expression experiment, 10M to 20M reads per
sample is usually enough.  If similar data exists it can be helpful to check
the read counts for key genes of interest to estimate the required depth.

Biological variability is usually the largest effect limiting the power of 
RNA-seq analysis.  The most improvement in an experiment will usually be 
achieved by increasing the biological replication to improve estimation of 
the biological variation.

It is often possible to design experiments where the analysis is done 
incrementally such that a pilot study is added to with an additional block of
samples or a pool of libraries is sequenced to additional depth. In these cases
care must be taken to balance the design in a manner that each stage is a
valid experiment in its own right.  This can allow a focused question to be 
answered in the first stage, with an ability to either address issues or 
progress to a second stage with additional questions.


## Sequencing options to consider

**How much total RNA is needed:**  
Many sequencing centres such as
[AGRF](http://www.agrf.org.au/resources/next-gen-resources/ngs-faqs#How much RNA do you need?)
recommend at least 250ng of total RNA for RNA sequencing. It is possible to go
as low as 100ng of total RNA, but results are not guaranteed. The quality of
RNA is also important when making libraries. A RNA Integrity Number (RIN) is a
number from 1 (poor) to 10 (good) and can indicate how much degradation there
is in the sample. A poor score can lead to over representation at the 3' end
of the transcript and low yield. Samples with low RIN scores (below 8) are
not recommended for sequencing.  Care should also be taken to ensure RIN is
consistent between conditions to avoid confounding this technical effect with
the biological question.

**Choosing an enrichment method:**  
Ribosomal RNA makes up >95% of total cellular RNA, so a preparation for
RNA-seq must either enrich for mRNA using poly-A enrichment, or deplete rRNA.
Poly-A enrichment is recommended for most standard RNA-seq experiments, but
will not provide information about microRNAs and other non-coding RNA species.
In general, ribo-depleted RNA-seq data will contain more noise, however, the
protocol is recommended if you have poor or variable quality of RNA as the 3’
bias of poly-A enrichment will be more pronounced with increased RNA
degradation. The amount of RNA needed for each method differs. For Poly-A
enrichment a minimum of 100ng is needed while for ribo-depletion, a minimum of
200ng is recommended.

**Choosing read type:**  
For basic differential expression analysis RNA-seq experiments, single-end
sequencing is recommended to obtain gene transcript counts. In more advanced
experiments, paired-ends are useful for determining transcript structure and
discovering splice variants.

**Choosing strandedness:**  

With a non-directional (unstranded) protocol, there is no way to identify
whether a read originated from the coding strand or its reverse complement.
Non-directional protocols allow mapping of a read to a genomic location, but
not the direction in which the RNA was transcribed. They are therefore used to
count transcripts for known genes, and are recommended for basic RNA-seq
experiments. Directional protocols (stranded) preserve strand information and
are useful for novel transcript discovery.

**Multiplexing:**  
Multiplexing is an approach to sequence multiple samples in the same
sequencing lane. By sequencing all samples in the same lane, multiplexing can
also minimise bias from lane effects.

**Spike-in controls:**  
RNA-seq spike-in controls are a set of synthetic RNAs of known concentration
which act as negative or positive controls. These controls have been used for
normalisation and quality control, but recent work has shown that the amount
of technical variability in their use dramatically reduces their utility.


## Summary

 - A good experimental design is vital for a successful experiment. If you're
   planning to work with a data analyst or bioinformatician, include them in
   the design stage.
 - Aim to minimise variability by identifying possible sources of variance in
   your samples.
 - Biological replicates are important. The number of biological replicates
   you should have should never be below 3. Technical replicates are often
   unnecessary.
 - Pilot studies are highly recommended for identifying how many replicates
   and how many reads you should have for enough statistical power in your
   experiment.
 - For basic RNA-seq experiments, poly-A enriched, single-ended, unstranded
   sequencing at depths of 10M to 20M is probably what you want.
