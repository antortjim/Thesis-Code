
# Index of what the scripts do

## PIPELINE SCRIPTS


## XIC_norm.py

Replicates the fraction normalization algorithm implemented in MaxQuant that aggregates XICs for the same peptide
collected across different fractions of the same sample.

## LFQ.py

Replicates the LFQ algorithm implemented in MaxQuant that estimates protein abundances by
minimising the global protein intensity differences based on the protein ratios.


### Input 

A `protein_ratios.txt` file, which is a matrix where

- Every row represents a protein group
- Every column represents a comparison between 2 samples
- The first column (no colname) shows the ids in the protein group 

This file can be obtained with process_peptide_intensities.R

### Output

A `protein_intensities.tsv` file, a table with one column per sample and one row per protein group
Every i,j cell represents the quantification of group i in sample j.



## process_peptide_intensities.R

Reads 


## moff_to_msqrob.R

Reads moff output in a peptide_summary_intensity_moFF_run.tab and the experimental design in `exp_name/data/experimental_design.tsv`.
Returns a data frame stored in the `exp_name/quantification/RSqM_signif.tsv` with the following fields:

- estimate (robust estimate of the log2FC of a contrast)
- se (standard error of the estimate)
- df (degrees of freedom)
- Tval (value of the student's t statistic)
- pval (p-value of the test)
- qval (corrected pvalue)
- signif (true if qval is less than 0.05)
- Protein.IDs (protein ids of the group)

This table can be fed to proteomic_analysis.R to create a volcano plot and perform gene set enrichment analysis

### Call
`nohup Rscript scripts/R/moff_to_msqrob.R --exp_name maxlfq --moff_file peptide_summary_intensity_moFF_run.tab --sample_filter "" --experiment_contrasts conditionH-conditionL --save_model &`

## proteomic_analysis.R

Reads MSqRob output from exp_name/quantification/RSqM_signif.tsv and performs differential abundance analysis, gsea and creates a volcano plot.

### Call
``

## 

## BENCHMARK SCRIPTS

## 