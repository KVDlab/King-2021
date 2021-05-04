Figure 7


Figure 7 is a visualization of TLR9 motifs.

1) We use python to perform most of the calculations. Including scrambling the sequences 1000x to generate a random distribution.
2) R is used for visulaization




1) Run 'dinucs.scrampled.py'. 
	This script uses '../shared_data/PaVE.complete-genomes.fas as input. It leverages 
		the functions in 'altschulEriksonDinuclShuffle.py' and compseq in Emboss to generate:
			1) ratio_sample-4.csv (contains ratios based on the actual sequences)
			2) ratio_test-4.csv (contains ratios based on 1000 randomly shuffled sequences)
			
3) Run Figure 7.R
	uses:
		../shared_data/4_PaVE.complete-genomes.fas_ratio.csv
		ratio_test-4.csv
		ratio_sample-4.csv
		
		
	generates:
		yang_4OE.csv (.csv containing the averaged Yang PV OE (observed/expected) values for each N-CG-N tetramer)
		4mer-comp_plots.csv (.csv containing the averaged Yang, Yin, and other OE (observed/expected) ratios for each N-CG-N tetramer)
		figure panels






