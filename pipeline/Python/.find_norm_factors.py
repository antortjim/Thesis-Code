# Given a tsv file containing protein ratios, this script performs LSQ optimization to reach protein intensities across samples
# Every row i of the protein ratios should correspond to the ith protein group
# Every column j should be the ratio of the jth combination

import pandas as pd
import numpy as np
from scipy.optimize import root
from tqdm import tqdm
import argparse
import os.path
import itertools
import subprocess
parser = argparse.ArgumentParser()
parser.add_argument("--root_dir", required=True)
parser.add_argument("--exp_name", required=True)
parser.add_argument("--n_peps", required=False, help="Set number of peptides to use in the normalization for debugging purposes", type=int)
arguments = vars(parser.parse_args())

output_moff_path = os.path.join(arguments["root_dir"], arguments["exp_name"], "peptideShaker_out", "PSM_reports", "output_moff_RAW")
output_moff = pd.read_csv(os.path.join(output_moff_path, "peptide_summary_intensity_moFF_run.tab"), sep = "\t")
experimental_design = pd.read_csv(os.path.join(arguments["root_dir"], arguments["exp_name"], "data", "experimental_design.tsv"), sep = "\t")
experimental_design.sort_values(by = ["Group", "Replicate"], inplace=True)

# A list of experiment + replicate combinations (H1, L1, H2, L2, etc)
exp_rep = experimental_design["Experiment"].map(str) + "_" + experimental_design["Replicate"].map(str)
# Dictionary storing the indices of each experiment_replicate in the experimental design
# i.e exp_rep_dic[H_1] is a list whith indices leading to the rows that belong to the fractions of H_1 (er_indices)
exp_rep_row = {e_r: np.where(exp_rep == e_r)[0] for e_r in exp_rep}

# Make it a list so that the len can be computed a priori
exp_rep_combinations = list(itertools.combinations(np.unique(exp_rep), 2))

n = output_moff.shape[0]
print("Found {} peptides".format(n))
if arguments["n_peps"] is not None:
    n = arguments["n_peps"]

print("Minimising for {} peptides".format(n))

xics_values = output_moff.iloc[:n, 2:].values

# def group_in_experiment(g, e):
#     exper = experimental_design["Experiment"].loc[experimental_design["Group"] == g]
#     return exper == e

def fun_xic_norm(N):
    # XICs for peptide P (all samples and fractions are in the same row)
    # I_{P,A} will be the sum of the product of XIC*N where condition = A
    # I_{P,A} is stored under key A of the I_P dictionary
    # H_P is the sum of the square of the absolute value of the logarithm of the ratio between all possible combinations of conditions
    # The combination is taken into account if neither I_PA nor I_PB are 0, otherwise it is ignored
    normalized_xics = xics_values * N

    # Create a dictionary I with a key for each experiment_replicate
    # The value of the key is set to the result of extracting the xic for samples belonging to e_r
    # and making a colSum. (np.sum(x, axis=1))
    # The value is then an array of length equal to the number of peptides where each ith cell
    # is the sum of the normalized intensities for peptide i in e_r
    I = {e_r: np.sum(normalized_xics[:, er_indices], axis=1) for e_r, er_indices in exp_rep_row.items()}

    # For every combination, take the logarithm of the intensity ratio for all peptides
    # Take the absolute value, square it and sum all the values.
    # This will be H for all peptides in the comparison c[0] / c[1].
    # We get one H for each comparison.
    H = np.sum(np.array([np.square(np.abs(np.log(I[c[0]] / I[c[1]]))) for c in exp_rep_combinations]).T, axis = 1)
    H = H[~np.isnan(H)]
    H = H[np.isfinite(H)]
    return H

#https://stackoverflow.com/questions/43995862/python-multivariate-non-linear-solver-with-constraints?utm_medium=organic&utm_source=google_rich_qa&utm_campaign=google_rich_qa
def constrainedFunction(x, f, lower, upper, minIncr=0.001):
     x = np.asarray(x)
     lower = np.asarray(lower)
     upper = np.asarray(upper)
     xBorder = np.where(x<lower, lower, x)
     xBorder = np.where(x>upper, upper, xBorder)
     fBorder = f(xBorder)
     distFromBorder = (np.sum(np.where(x<lower, lower-x, 0.)) + np.sum(np.where(x>upper, x-upper, 0.)))
     return (fBorder + (fBorder + np.where(fBorder>0, minIncr, -minIncr))*distFromBorder)

def main(minimize=True):
    n_samples = xics_values.shape[1]
    x0 = np.full(n_samples, init_value)
    if minimize:
        # find the normalization fators that minimise fun_xic_norm
        N = root(constrainedFunction, x0=x0, method="lm", args=(fun_xic_norm, np.array([0,]*n_samples), np.array([20,]*n_samples))).x
        print("Saving norm factors to file")
        np.savetxt(fname=os.path.join(output_moff_path, "normalization_factors.txt"), X=N, fmt="%10.5f", delimiter="\t")


    else:
        N = x0

if __name__ == "__main__":
    init_value = 1
    main(True)
    main(False)
