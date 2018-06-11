# Given a tsv file containing protein ratios, this script performs LSQ optimization to reach protein intensities across samples
# Every row i of the protein ratios should correspond to the ith protein group
# Every column j should be the ratio of the jth combination

import pandas as pd
import numpy as np
from scipy.optimize import least_squares
from tqdm import tqdm
import argparse
import os.path

parser = argparse.ArgumentParser()
parser.add_argument("--input_dir", required=True)
parser.add_argument("--output_dir", required=True)
arguments = vars(parser.parse_args())


protein_ratios = pd.read_csv(os.path.join(arguments["input_dir"], "protein_ratios.tsv"), sep = "\t").iloc[:, :15]
combinations = np.array(list(map(lambda x: x.split("/"), protein_ratios.columns.values.tolist())))
n = protein_ratios.shape[0]
sample_names = np.sort(np.unique(combinations.flatten()))
mapper = {k: i for i, k in enumerate(sample_names)}

def map_pos(c):
    rm = mapper
    pos = [rm[cc] for cc in c]
    return(pos)

def fun_lfq(x, i):
    ratio = protein_ratios.values[i,:]
    
    valid_ratios = ratio != 0
    c1 = combinations[:,0]
    c2 = combinations[:,1]
#    print(x[map_pos(c1)])
#    print(x[map_pos(c2)])
#    print(valid_ratios)
    to_sum = (np.log(ratio) - np.log(x[map_pos(c1)]) + np.log(x[map_pos(c2)]))
    to_sum = to_sum[valid_ratios]
    return np.sum(to_sum)


if __name__ == "__main__":
    init_value = 1
    protein_intensities = np.zeros((n, len(sample_names)))
    for i in tqdm(range(n)):
        minimisation = least_squares(fun_lfq, x0 = np.array([init_value] * len(sample_names)), kwargs={"i": i}).x
        protein_intensities[i,:] = np.where(minimisation == init_value, 0, minimisation)
    

    np.savetxt(fname = os.path.join(arguments["output_dir"], "protein_intensities.tsv"), X=protein_intensities, fmt = "%10.5f", delimiter="\t", header="\t".join(sample_names.tolist()))
