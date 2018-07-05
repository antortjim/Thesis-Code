# Pipeline documentation #

 * [Introduction](#introduction)
 * [Minimum Requirements](#minimum-requirements)
 * [Scripts](#scripts)
 * [Output Data](#output-data)
---

## Introduction ##

This repo provides a complete open source, free and (almost) 100 % Linux supported pipeline for the label-free quantification analysis of RAW MS spectra files.

It is achieved by performing the main following tasks:

- Search spectra against protein database and find Peptide-to-Spectrum-Matches (PSMs) using one or several search engines.
- Validation and FDR filtering of the search results.
- Match between runs and apex intensity extraction.
- Quantification proper.

## Minimum requirements ##

## Scripts ##


 * [Load settings](#create_settings_file.sh)
 * [Create decoy database](#create_decoy_database.sh)
 * [Search spectra](#search_all_mgf.sh)
 * [Integrate and validate results](#call_peptide_shaker.sh)
 * [MBR and Apex intensity extraction](#call_moFF.sh)
 * [Relative quantification](#moff_to_msqrob.R)
 * [Biological inference](#proteomic_analysis.R)
 

 
## create_settings_file.sh

## create_decoy_database.sh

## search_all_mgf.sh

## call_peptide_shaker.sh

## call_moFF.sh


## moff_to_msqrob.R

Reads moFF-like output in a peptide_summary_intensity_moFF_run[_SUFFIX].tab and the experimental design in `exp_name/data/experimental_design.tsv`.
Returns a robust estimation of the log2(FC) for as many proteins as possible, together with test statistics

This table can be fed to proteomic_analysis.R to create a volcano plot and perform gene set enrichment analysis

If the `--fraction_normalized` flag is passed, the moFF file is assumed to be produced by `aggregate_fractions` and not moFF itself.
This means that the experimental annotation needs to be edited to remove the fractions
The value of the `--suffix` is a string appended to the RSqM_signif filename
If `--save_model` the session is saved to the export folder, set by default to the quantification folder of the experiment

#### Details

The script makes use of the MSqRob package for robust protein quantification of MaxQuant/moFF output.
In turn this package makes use of the data structures implemented in the MSnBase package (MSnSet data).

It works by:

1. Importing the data with `import2MSnSet`.
2. Preprocess with `preprocess_MSnSet`.
3. Compile to a protdata object (implemented in MSqRob).
4. Fit a ridge regression model with Huber weights and empirical bayes estimation of the variance (robust regression).
5. Build a contrast matrix and perform hypothesis testing with `test.protLMcontrast`.
6. Adjust results for multiple testing with `prot.p.adjust`.


### Input

* A `peptide_summary_intensity_moFF_run.tab` file created by moFF or by aggregate_fractions.R
* A `experimental_design.tsv` file storing the sample organisation
* Several flags to fine tune the behaviour: which contrast to test, etc.

### Output

A `RSqM_signif` file saved to the default export folder. It is a table with the following fields:

- estimate (robust estimate of the log2FC of a contrast)
- se (standard error of the estimate)
- df (degrees of freedom)
- Tval (value of the student's t statistic)
- pval (p-value of the test)
- qval (corrected pvalue)
- signif (true if qval is less than 0.05)
- Protein.IDs (protein ids of the group)

### Call

`nohup Rscript scripts/R/moff_to_msqrob.R --root_dir `pwd` --exp_name maxlfq --moff_file peptide_summary_intensity_moFF_run.tab --sample_filter "" --experiment_contrasts conditionH-conditionL --save_model --suffix "" [--fraction_normalized] &`

## proteomic_analysis.R

Reads MSqRob output from `exp_name/quantification/RSqM_signif.tsv` and performs differential abundance analysis, Protein Set Enrichment Analysis (PSEA, similar to GSEA) and creates a volcano plot.

## Input

* RSqM_signif.tsv file created by moff_to_msqrob.R 

## Output

Plots:

-Volcano plot showing for every protein group:
  - on the x axis the value of the estimated log2FC
  - on the y axis the value of the -log10 p-value

-Histogram showing the estimate distribution

-Lists of up-regulated and down-regulated proteins.
-List of proteins with significant change but beneath the log2fc_threshold.
-List of proteins for which the estimate could not be computed. This is due to most peptides being missing in all replicates of one of the conditions of the contrast
  
Plots and summary statistics


### Call

`nohup Rscript scripts/R/proteomic_analysis.R --root_dir `pwd` --exp_name maxlfq --suffix "" --log2fc_threshold 1 &`


## Peptide and protein search and inference


## find_normalization_factors.py

Replicates the "fraction normalization" algorithm implemented in MaxQuant that aggregates XICs for the same peptide
collected across different fractions of the same sample. This is done by performing Levenberg-Marquandt (LM) minimisation
of the squared logarithm of the peptide intensity ratio across any two samples for all samples and peptides

### Input

* A `peptide_summary_intensity_moFF_run.tab` file created by moFF
* A `experimental_design.tsv` file storing the sample organisation

It takes an optional argument `--n_peps` if the user wants to minimise on less peptides than those available for debugging / speed up purposes.

### Output

A normalisation_factors.txt file storing the normalisation factor for sample ith in the ith row

## Call

`python scripts/Python/find_norm_factors.py --root_dir `pwd` --exp_name maxlfq --n_peps 1000`

## aggregate_fractions.py

Completes the normalization process initiated by `find_norm_factors.py` by taking the estimated root of the LM minimisation (the normalization factors)
and aggregates  in order to remove the fraction effect.

### Input

* A `peptide_summary_intensity_moFF_run.tab` file created by moFF
* A `experimental_design.tsv` file storing the sample organisation
* A `normalization_factors.txt` file storing the solution to the LM minimisation

### Output

* A `peptide_summary_intensity_moFF_run_fraction_normalized.tab` file with format identical to that created by moFF but now featuring one column per experiment+replicate.
The name of the column is built like this "{}{}".format(experiment, replicate). So the column representing data from all fractions of the sample from condition L and replicate 1
will be called L1.

## Call

`python scripts/Python/aggregate_fractions.py --root_dir `pwd` --exp_name maxlfq`


### PEPTIDE AGGREGATION

## LFQ.py

Replicates the LFQ algorithm implemented in MaxQuant that estimates protein abundances by minimising the global protein intensity differences based on the protein ratios.


### Input 

A `protein_ratios.txt` file, which is a matrix where

- Every row represents a protein group.
- Every column represents a comparison between 2 samples.
- The first column (no colname) shows the ids in the protein group.

This file can be obtained with process_peptide_intensities.R

### Output

A `protein_intensities.tsv` file, a table with one column per sample and one row per protein group
Every i,j cell represents the quantification of group i in sample j.





## BENCHMARK SCRIPTS

## 