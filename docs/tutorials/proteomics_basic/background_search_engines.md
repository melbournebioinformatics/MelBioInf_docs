#How Search Engines Work

When choosing search engine parameters it can be helpful to understand the basic algorithms that most search engines share. The workflow of a typical search engine is roughly as outlined below.

1. The next spectrum in the dataset is loaded for consideration
2. The spectrum parent mass and the parent mass tolerance is used to select a small set of matching peptides from the database.  This is a crucial step because the number of peptides that fall within the matching mass window will determine the effective size of the search space.  Database search space size also depends on many other factors including
	- The size of the protein database
	- Variable modifications allowed on peptides
	- Parent ion mass tolerance
	- Number of allowed missed enzymatic cleavages
	- Specificity of the enzyme used for digestion
3. Each of the matching peptides is scored against the spectrum. The nature of this scoring is where search engines typically differ from each other.
4. The highest scoring `peptide spectrum match` (PSM) is recorded along with its score.
5. Some form of global analysis of (PSM) scores is performed in order to determine a threshold of significance

There are many excellent presentations online that explain this in more detail.  Although it's old, I recommend [this presentation by Brian Searle](http://www.slideshare.net/ProteomeSoftware/interpreting-msms-results)

