import argparse
import pandas as pd

parser = argparse.ArgumentParser()
parser.add_argument("--group", required=True)
parser.add_argument("--exp_design", required=True)
parser.add_argument("--prepend", required=True)
parser.add_argument("--append", required=True)
parser.add_argument("--ps_out")
args = parser.parse_args()
arguments = vars(args)

exp_design = pd.read_csv(arguments["exp_design"], sep="\t", header=0)
names = exp_design.loc[exp_design["Group"] == arguments["group"],:]["Name"].values
sample_names = [arguments["prepend"] + "/" + e + arguments["append"] for e in names]
output = " ".join(sample_names) 
print(output)

