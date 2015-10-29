#Protein Databases

In a perfect experiment we would obtain fragment ions for all the `b,y` pairs of each peptide.  If peaks can be unambiguously identified for all these pairs then the sequence of a peptide can simply be read off from the fragmentation spectrum itself.  Unfortunately this is almost never the case using current instrumentation, and the only practical method to determine the sequences of peptides and proteins present in a sample is to compare spectra with a database of potential proteins.  This `database` is usually just a `FASTA` formatted file containing amino acid sequences for all known proteins from your study organism.  Constructing this database is a crucial first step in any proteomics analysis because only peptide sequences present in the database will appear in the results.  In order to detect a peptide, its exact sequence must be explicitly included in the database. 

## Large vs Small Database

Since it is impossible to detect a peptide unless it is present in the search database, one might consider using a very large database such as the full content of `NCBInr`.  There are two problems with this.  The most important is that the sensivity of a search goes down as the search space goes up, which means that searches on large databases often return far fewer hits.  Another, more practical issue is that searching a large database often takes an extremely long time, and might even crash your search engine. 

Note that very small databases can also cause problems.  In particular, some search engines, and most search engine post-processing statistical tools attempt to model the shape of peptide-sequence-match (PSM) score distributions.  With a very small database (or with very few spectra) it may not be possible to model these distributions accurately.  In most practical situations this is not an issue.

## Typical sources of data for search databases

**Uniprot.org**: This is the canonical resource for publicly available protein sequences. It includes two large databases `SwissProt`, which contains manually curated sequences and `Trembl` which contains sequences automatically generated from genomic and transcriptomic data.

**An Organism Specific database**: In some cases, a community of researchers working on specific organisms will create their own sequence data repositories.  Some of these are well maintained and are the best source of data for that study organism.  Examples include `PlasmoDB` for Malaria for `Flybase` for drosophila.

**Transcriptome derived sequences**: If you are working on an organism for which public sequence data are scarce, it may be worth obtaining transcriptomic sequences for it.  If sufficient data are obtained, the resulting assembled transcript sequences can be translated to form a good quality proteomic database.

**Other**: Depending on the project you might want to include sequences for specific variants of interest, or a six-frame translation of a genome. 

## Should I include decoys?

Decoys are often useful, but not always needed.  Often, the decision to include decoys depends on the requirements of software that is used downstream of the search. Examples on this wiki that make use of `Peptide Prophet` typically use decoys because it can use these to 'pin down' the negative distribution.

