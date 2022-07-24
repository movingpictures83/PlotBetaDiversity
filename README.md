# PlotBetaDiversity
# Language: R
# Input: TXT
# Output: PREFIX
# Tested with: PluMA 1.1, R 4.0.0
# Dependencies: dplyr 1.0.8, microbiome 1.12, vegan 2.5.7, ggpubr 0.4.0, Matrix 1.3.2, reshape2 1.4.4, ggplot2 3.3.5

PluMA plugin that computes Beta Diversity plots, one for pairwise distance within groups, 
and the other for differential analysis.

The following are specified in the input TXT file, as tab-delimited keyword-value pairs.

otufile: OTU abundances (CSV)
mapping: Mapping table (CSV)
tree: Phylogenetic tree (CSV)
distance: Metric for distance calculation (STRING)
differential: Metric for differential analysis (STRING).
column: The column to differentiate
colors: The colors for differentiation

The plot will then be generated using the provided prefix (prefix.pdf), with a supplemental prefix.csv file containing diversity values.
