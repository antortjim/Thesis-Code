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
arguments = vars(parser.parse_args())


output_moff_path = os.path.join(arguments["root_dir"], arguments["exp_name"], "peptideShaker_out", "PSM_reports", "output_moff_RAW")

### OCCAM RAZOR here using the R function from MSqRob
#occam_path = os.path.join(arguments["root_dir"], "scripts", "bash", "occam_razor.sh")
#subprocess.call("{} {} {} {}".format(occam_path, arguments["root_dir"], arguments["exp_name"], output_moff_path), shell=True)
###########################################
output_moff = pd.read_csv(os.path.join(output_moff_path, "peptide_summary_intensity_moFF_run.tab"), sep = "\t")


experimental_design = pd.read_csv(os.path.join(arguments["root_dir"], arguments["exp_name"], "data", "experimental_design.tsv"), sep = "\t")
experimental_design.sort_values(by = ["Group", "Replicate"], inplace=True)

groups = experimental_design["Group"].values
experiment = experimental_design["Experiment"].values
replicate = experimental_design["Replicate"].values
groups_unique = np.unique(groups)

# A list of experiment + replicate combinations (H1, L1, H2, L2, etc)
exp_rep = experimental_design["Experiment"].map(str) + "_" + experimental_design["Replicate"].map(str)
# Dictionary storing the indices of each experiment_replicate in the experimental design
# i.e exp_rep_dic[H_1] is a list whith indices leading to the rows that belong to the fractions of H_1
exp_rep_dict = {e_r: np.where(exp_rep == e_r)[0] for e_r in exp_rep}
group_indices = {g: np.where(g == groups)[0].tolist() for g in groups_unique}

for key, value in group_indices.items():
    new_dict = {}
    for v in value:
        r = experimental_design["Replicate"].iloc[v]
        new_dict[r] = v
    group_indices[key] = new_dict

# Make it a list so that the len can be computed a priori
exp_rep_combinations = list(itertools.combinations(np.unique(exp_rep), 2)) 

n = output_moff.shape[0]
xics = output_moff.iloc[:n, 2:]

def group_in_experiment(g, e):
    exper = experimental_design["Experiment"].loc[experimental_design["Group"] == g]
    return exper == e

def fun_xic_norm(N):
    # XICs for peptide P (all samples and fractions are in the same row)
    # I_{P,A} will be the sum of the product of XIC*N where condition = A
    # I_{P,A} is stored under key A of the I_P dictionary
    # H_P is the sum of the square of the absolute value of the logarithm of the ratio between all possible combinations of conditions
    # The combination is taken into account if neither I_PA nor I_PB are 0, otherwise it is ignored
    normalized_xics = xic_values * N
    
    I = {e_r: np.sum(normalized_xics[:, indices], axis=1) for e_r, indices in exp_rep_dict.items()}
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
    x0 = np.full(xics.shape[1], init_value)
    if minimize:
        # find the normalization fators that minimise fun_xic_norm
        N = root(constrainedFunction, x0=x0, method="lm", args=(fun_xic_norm, np.array([0,]*xics.shape[1]), np.array([20,]*xics.shape[1]))).x
        prefix = "_fraction_normalization"
        print("Saving norm factors to file")
        np.savetxt(fname=os.path.join(output_moff_path, "normalization_factors.tsv"), X=N, fmt="%10.5f", delimiter="\t")


    else:
        prefix = "_raw_sum"
        N = x0

    intensities = xics.values * N 

    I = {}
    result = []
    k = 0
    colnames = []
    for r in np.unique(replicate):
        for e in np.unique(experiment):
             print(e)
             II = []
             for i, g in enumerate(groups_unique):
                 if experiment[i] == e and r in group_indices[g].keys():
                     new_key = "{}_{}".format(g, r)
                     col = group_indices[g][r]
                     inten = intensities[:, col]
                     I[new_key] = inten 
                     II.append(inten)
             new_key = e + "_" + str(r)
             II = np.vstack(II).T
             II = np.sum(II, axis=1)
             colnames.append("{}_{}".format(e, r))
             result.append(II)
                     
    result = np.vstack(result).T
    result = pd.DataFrame.from_dict(result, orient="columns")
    result.columns = colnames
    result = pd.concat([xics_full.iloc[:, :2], result], axis = 1)
    result.columns = ["Sequence", "Protein.IDs"] + colnames
    print("Saving aggregated intensities to file")
    result.to_csv(os.path.join(arguments["root_dir"], arguments["exp_name"], "peptideShaker_out", "PSM_reports", "output_moff_RAW", "peptide_summary_intensity_moFF_run" + prefix + ".tab"),
                  sep="\t", index=False)

if __name__ == "__main__":
    init_value = 1
    xic_values = xics.values
 
    main(True)
    main(False)
