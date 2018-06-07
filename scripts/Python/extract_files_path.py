import argparse
import pandas as pd

parser = argparse.ArgumentParser()
parser.add_argument("--sample", required=True)
parser.add_argument("--exp_design", required=True)
parser.add_argument("--prepend", required=True)
parser.add_argument("--append", required=True)
parser.add_argument("--ps_out")
args = parser.parse_args()
arguments = vars(args)

exp_design = pd.read_csv(arguments["exp_design"], sep="\t", header=0)
files = exp_design.loc[exp_design["Experiment"]== arguments["sample"],:]["Name"].values
sample_names = [arguments["prepend"] + "/" + e.replace(".mgf", "") + arguments["append"] for e in files]
output = " ".join(sample_names) 
print(output)

