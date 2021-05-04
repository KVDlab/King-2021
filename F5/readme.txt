Figure 5


Figure 5 is a visualization of codon usage.

1) We use python to remove overlapping parts of the viral ORFs from the PV genomes. 
	These edited ORFs are assembled into the 'concatenated.fas' file.
2) We calculate the codon usage metrics in figure 5A
3) We also calculate the amino acid usage based on each of these sequences 
4) R calculates the RSCU values and generates the figure panels




1) Run 'delete_overlaps.py'. 
	This script uses '../shared_data/reference.csv' and ../shared_data/PaVE.csv to generate:
		E6.fas, E7.fas, E1.fas, E2.fas, L1.fas, L2.fas, and concatenated.fas.
		The file concatenated.fas will be used for the next steps

2) Run 'cusp.py'. 
	This script leverages cusp and codcmp from Emboss to generate and compare codon usage tables
		final output: 'input4R.csv'
	We also calculate amino acid usage as described https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2590925/
		output: 'aa_composition.csv' and 'aa_composition_values.csv'

3) Run Figure 5.R
	generates 'rscu.csv', 'rscu_plots.csv' and the individual figure panels






