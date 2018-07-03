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
parser.add_argument("--suffix", required=True)
parser.add_argument('--norm_factors', dest='norm_factors', action='store_true')
parser.add_argument('--no_norm_factors', dest='norm_factors', action='store_false')
parser.set_defaults(feature=True)
arguments = vars(parser.parse_args())


output_moff_path = os.path.join(arguments["root_dir"], arguments["exp_name"], "peptideShaker_out", "PSM_reports", "output_moff_RAW")
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

# dictionary storing for each group, a new dict with keys set to the available
# replicates for that group. The value under each of these keys will be the
# index of the corresponding row in experimental_design
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
xics_values = output_moff.iloc[:n, 2:].values
n_samples = xics_values.shape[1]

if arguments["norm_factors"]:
    N = np.loadtxt(fname=os.path.join(output_moff_path, "normalization_factors.txt"), delimiter="\t")
else:
    N = np.full(n_samples, 1)

def aggregate_fractions(xics_values, N):
    # Multiply each intensity with the normalization factor
    intensities = xics_values * N

    # Now sum the different fractions
    #I = {}
    result = []
    k = 0
    colnames = []

    # For each replicate in each experiment
    for r in np.unique(replicate):
        for e in np.unique(experiment):

             # Initialize a list storing the normalized intensities for the current experiment_replicate (sample)
             II = []
             for i, g in enumerate(groups_unique):
                 # check that the ith experiment is the current experiment
                 # and that the current experiment has replicate r
                 if experiment[i] == e and r in group_indices[g].keys():
                     # the column we neet to extract from the intensities table
                     # is equal to the index of the rth replicate in the g group
                     col = group_indices[g][r]
                     inten = intensities[:, col]
                     #I[new_key] = inten
                     II.append(inten)
             # once the normalized intensity of all the fractions
             # under experiment e and replicate r
             # are stored in a list, make the
             II = np.vstack(II).T
             II = np.sum(II, axis=1)
             colnames.append("sumIntensity_{}{}".format(e, r))
             result.append(II)

    result = np.vstack(result).T
    result = pd.DataFrame.from_dict(result, orient="columns")
    result.columns = colnames
    result = pd.concat([output_moff.iloc[:, :2], result], axis = 1)
    result.columns = ["peptide", "prot"] + colnames
    print("Saving aggregated intensities to file")
    peptide_summary_file = "peptide_summary_intensity_moFF_run_fraction_normalized{}.tab".format(arguments["suffix"])

    result.to_csv(os.path.join(arguments["root_dir"], arguments["exp_name"], "peptideShaker_out", "PSM_reports", "output_moff_RAW", peptide_summary_file),
                  sep="\t", index=False)

aggregate_fractions(xics_values, N)
