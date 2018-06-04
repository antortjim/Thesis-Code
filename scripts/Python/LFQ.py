# Given a tsv file containing protein ratios, this script performs LSQ optimization to reach protein intensities across samples
# Every row i of the protein ratios should correspond to the ith protein group
# Every column j should be the ratio of the jth combination

import pandas as pd
import numpy as np
from scipy.optimize import least_squares
from tqdm import tqdm
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--input_dir", required=True)
parser.add_argument("--output_dir", required=True)
arguments = vars(parser.parse_args())


protein_ratios = pd.read_csv(os.path.join(arguments["input_dir"], "protein_ratios.tsv"), sep = "\t").iloc[:,:15]
combinations = np.array(list(map(lambda x: x.split("/"), protein_ratios.columns.values.tolist())))
keys = np.unique(combinations).tolist()
n = protein_ratios.shape[0]


def map_pos(c):
    rm = {"H1": 0, "H2": 1, "H3": 2, "L1" : 3, "L2": 4, "L3": 5}
    return [rm[cc] for cc in c]

def fun_lfq(x, i):
    ratio=protein_ratios.values[i,:]
    
    valid_ratios = ratio != 0
    c1 = combinations[:,0]
    c2 = combinations[:,1]
    to_sum = (np.log(ratio) - np.log(x[map_pos(c1)]) + np.log(x[map_pos(c2)]))[valid_ratios]
    return np.sum(to_sum)


if __name__ == "__main__":
    init_value = 1
    protein_intensities = np.zeros((n, 6))
    for i in tqdm(range(n)):
        minimisation = least_squares(fun_lfq, x0 = np.array([init_value]*6), kwargs={"i": i}).x
        protein_intensities[i,:] = np.where(minimisation == init_value, 0, minimisation)
    

    np.savetxt(fname = os.path.join(output_dir, "protein_intensities.tsv"), X=protein_intensities, fmt = "%10.5f", delimiter="\t")
