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
files = exp_design.loc[exp_design["sample"] == arguments["sample"],:]["file"].values

sample_names = [arguments["prepend"] + "/" + e.replace(".mgf", "") + arguments["append"] for e in files]


if arguments["append"] == "_Default_PSM_Report.txt":
    sample_mapping = pd.read_csv(arguments["prepend"] + "/../" + "sample_mapping.txt")
    ps_sample_names = []
    for sn in sample_names:
        ps_sn = sample_mapping.loc[sample_mapping.iloc[:,0] == sn,:].iloc[:,1]
        ps_sample_names.append(ps_sn)
  
    sample_names = ps_sample_names


output = " ".join(sample_names) 
print(output)

