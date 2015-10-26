# De novo genome assembly using Velvet

## Background

### Introduction to de novo assembly

DNA sequence assembly from short fragments (< 200 bp) is often the first step of any bioinformatic analysis. The goal of assembly is to take the millions of short reads produced by sequencing instruments and re-construct the DNA from which the reads originated.

The sequence assembly issue was neatly summed up by the following quote:

*"The problem of sequence assembly can be compared to taking many copies of a book, passing them all through a shredder, and piecing a copy of the book back together from only shredded pieces.  The book may have many repeated paragraphs,  and some shreds may be modified to have typos.  Excerpts from another book may be added in, and some shreds may be completely unrecognizable."* – Wikipedia: Sequence assembly.

An addition to the above for paired end sequencing is that now some of the shreds are quite long but only about 10% of the words from both ends of the shred are known.

This tutorial describes de novo assembly of Illumina short reads using the Velvet assembler (Zerbino  et al. 2008, 2009) and the Velvet Optimiser (Gladman & Seemann, 2009) from within the Galaxy workflow management system.

----------------------------------

## The Galaxy workflow platform

### What is Galaxy?

Galaxy is an online bioinformatics workflow management system. Essentially, you upload your files, create various analysis pipelines and run them, then visualise your results.

Galaxy is really an interface to the various tools that do the data processing; each of these tools could be run from the command line, outside of Galaxy. Galaxy makes it easier to link up the tools together and visualise the entire analysis pipeline.

Galaxy uses the concept of 'histories'. Histories are sets of data and workflows that act on that data. The data for this workshop is available in a shared history, which you can import into your own Galaxy account

Learn more about Galaxy here

**Figure 1: The Galaxy interface**

Tools on the left, data in the middle, analysis workflow on the right.

<img src="../media/screenshot-fastqc-output.png" width=100% />

------------------------------

## De novo assembly with Velvet and the Velvet Optimiser.

### Velvet

Velvet is software to perform dna assembly from short reads by manipulating de Bruijn graphs. It is capable of forming long contigs (n50 of in excess of 150kb) from paired end short reads. It has several input parameters for controlling the structure of the de Bruijn graph and these must be set optimally to get the best assembly possible. Velvet can read Fasta, FastQ, sam or bam files. However, it ignores any quality scores and simply relies on sequencing depth to resolve errors. The Velvet Optimiser software performs many Velvet assemblies with various parameter sets and searches for the optimal assembly automatically.

### de Bruijn graphs

A de Bruijn graph is a directed graph which represents overlaps between sequences of symbols. The size of the sequence contained in the nodes of the graph is called the word-length or k-mer size. In Figure 2, the word length is 3. The two symbols are 1 and 0. Each node in the graph has the last two symbols of the previous node and 1 new symbol. Sequences of symbols can be produced by traversing the graph and adding the “new” symbol to the growing sequence.

**Figure 2: A de Bruijn graph of word length 3 for the symbols 1 and 0.**

<img src="../media/debruijn.jpg" width=80% />

*From: https://cameroncounts.wordpress.com/2015/02/28/1247/*

Velvet constructs a de Bruijn graph of the reads. It has 4 symbols (A, C, G and T - N’s are converted to A’s) The word length (or k-mer size) is one of Velvet’s prime parameters.

Velvet is not the only assembly software that works in this manner. Euler, Edena and SOAP de novo are examples of others.

### The Velvet algorithm

#### *Step 1: Hashing the reads.*

* Velvet breaks up each read into k-mers of length k.
  * A k-mer is a k length subsequence of the read.
  * A 36 base pair long read would have 6 different 31-mers.
* The k-mers and their reverse complements are added to a hash table to categorize them.
* Each k-mer is stored once but the number of times it appears is also recorded.
* This step is performed by “velveth” - one of the programs in the Velvet suite.

#### *Step 2: Constructing the de Bruijn graph.*

* Velvet adds the k-mers one-by-one to the graph.
* Adjacent k-mers overlap by k-1 nucleotides.
* A k-mer which has no k-1 overlaps with any k-mer already on the graph starts a new node.
* Each node stores the average number of times its k-mers appear in the hash table.
* Figure 3 shows a section of a de Bruijn graph constructed by Velvet for k=5.
* Different sequences can be read off the graph by following a different path through it. (Figure 3)

**Figure 3: Section of a simple de Bruijn graph of reads with k-mer size 5. Coloured sequences are constructed by following the appropriately coloured line through the graph.** *(Base figure Zerbino et al 2008.)*

<img src="../media/velvet-graph.png" width=80% />

#### *Step 3: Simplification of the graph.*

* Chain merging: When there are two connected nodes in the graph without a divergence, merge the two nodes.
* Tip clipping: Tips are short (typically) chains of nodes that are disconnected on one end. They will be clipped if their length is < 2 x k or their average k-mer depth is much less than the continuing path.
* Bubble removal: Bubbles are redundant paths that start and end at the same nodes (Figure 4.) They are created by sequencing errors, biological variants or slightly varying repeat sequences.
  * Velvet compares the paths using dynamic programming.
  * If they are highly similar, the paths are merged.
* Error removal: Erroneous connections are removed by using a “coverage cutoff”. Genuine short nodes which cannot be simplified should have a high coverage. An attempt is made to resolve repeats using the “expected coverage” of the graph nodes.
* Paired end read information: Velvet uses algorithms called “Pebble” and “Rock Band” (Zerbino et al 2009) to order the nodes with respect to one another in order to scaffold them into longer contigs.

**Figure 4: Representation of “bubbles” in a Velvet de Bruijn graph.**  *(Base figure Zerbino et al 2008.)*

<img src="../media/velvet-bubbles.png" width=80% />

#### *Step 4: Read off the contigs.*

* Follow the chains of nodes through the graph and “read off” the bases to create the contigs.
* Where there is an ambiguous divergence/convergence, stop the current contig and start a new one.

### K-mer size and coverage cutoff values

The size of the k-mers that construct the graph is very important and has a large effect on the outcome of the assembly. Generally, small k-mers create a graph with increased connectivity, more ambiguity (more divergences) and less clear “paths” through the graph. Large k-mers produce graphs with less connectivity but higher specificity. The paths through the graph are clearer but they are less connected and prone to breaking down.

The coverage cutoff c used during the error correction step of Velvet also has a significant effect on the output of the assembly process. If c is too low, the assembly will contain nodes of the graph that are the product of sequencing errors and misconnections. If c is too high, it can create mis-assemblies in the contigs and destroys lots of useful data.

Each dataset has its own optimum values for the k-mer size and the coverage cutoff used in the error removal step. Choosing them appropriately is one of the challenges faced by new users of the Velvet software.

### Velvet Optimiser

The Velvet Optimiser chooses the optimal values for k and c automatically by performing many runs of Velvet (partially in parallel) and interrogating the subsequent assemblies.  It uses different optimisation functions for k and c and these can be user controlled.

It requires the user to input a range of k values to search (to cut down on running time).

## References
http://en.wikipedia.org/wiki/Sequence_assembly

Zerbino DR, Birney E, Velvet: algorithms for de novo short read assembly using de Bruijn graphs, Genome Research, 2008, 18:821-829

Zerbino DR, McEwen GK, Margulies EH, Birney E, Pebble and rock band: heuristic resolution of repeats and scaffolding in the velvet short-read de novo assembler. PLoS One. 2009; 4(12):e8407.

Gladman SL, Seemann T, Velvet Optimiser, http://www.vicbioinformatics.com/software.shtml 2009.
