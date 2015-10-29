# Protein Prophet

The development of the Protein Prophet statistical models, and its associated program was a big step forward for practitioners wanting to perform automated, large scale protein inference.  The original paper describing protein prophet is worth reading. It's citation is;

>Nesvizhskii, A. I., Keller, A., Kolker, E. & Aebersold, R. A Statistical Model for Identifying Proteins by Tandem Mass Spectrometry. Anal. Chem. 75, 4646–4658 (2003).

The practical reality of using Protein Prophet is a little different however as the program has undergone several significant developments since its original publication, and interpretation of Protein Prophet groupings can be challenging.  Let's look at a few examples;

**Uniquely identified protein**

> sp|P00761|TRYP_PIG

Search for this protein in the Protein Prophet results file. It should have `group_probability` and `protein_probability` scores of `1.0`.  All of the three peptides that contribute evidence for this protein map uniquely to this protein alone so there are no other proteins in this group.

**Indistinguishable Protein**

>sp|O08600|NUCG_MOUSE

In this case there still just one entry for the protein, but `Protein Prophet` lists another protein `tr|Q3UN47|Q3UN47_MOUSE` in the indistinguishable proteins column.  This protein is indistinguishable from the primary entry `sp|O08600|NUCG_MOUSE` because all of the identified peptides are shared between both.

**A well behaved protein group**

>sp|O08677|KNG1_MOUSE

This protein is part of a smallish group of similar proteins.  The overall group probability is high (`1.0`) but probabilities group members are different.  The first member of the group has a high probability `0.99` but all other members have probabilities of `0.0`.  This is because all of the high scoring peptides are contained in the first entry.  Evidence for the other entries consists of either (a) peptides that are contained in the first entry or (b) peptides with very low scores. Protein Prophet uses the principle of Occam's razor;

>plurality should not be posited with out necessity

In other words, unless otherwise indicated by a unique peptide, we should assume that shared peptides come from the same protein.

**Anomalous groups**

In rare cases Protein Prophet fails produces strange results when its algorithm fails to converge.  This can result in situations where the group probability is high (1.0) but all of the member proteins within the group are assigned a probability of 0. 


