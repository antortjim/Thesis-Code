import pymc3 as pm
import theano
import pandas as pd
import os.path
import glob
import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm
np.random.seed(123)
import seaborn as sns
sns.set_style('whitegrid')
from BayesQuant import BayesQuant 
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--n', type=int, required=True,
                    help='number of proteins to analyze per peptide_n')
parser.add_argument('--top', type=int, required=True,
                    help='max number of peptides to use. If more are available, they are filtered')

parser.add_argument('--input', type=str, required=True,
                    help="path to m1_intensities file")

parser.add_argument('--output', type=str, required=True,
                    help="output file")

args = vars(parser.parse_args())

n=args["n"]
top=args["top"]

bayesquant = BayesQuant(plot_dir="thp1/plots/",traces_dir="thp1/traces")
data,_ = bayesquant.read_data(data_path=args["input"], features_path=None)
counts0 = np.unique(data.protein, return_counts=True)
protein_counts = {str(c): counts0[0][counts0[1] == c] for c in np.unique(counts0[1])}

if os.path.isfile(args["output"]):
    bq = pd.read_csv(args["output"], sep="\t", index_col=0)
else:
    bq = pd.DataFrame({"mean": [], "sd": [], "hpd_2.5":[], "hpd_97.5":[],"Organism":[], "n_peptides":[]})


for n_pep in range(2,6):
    proteins = protein_counts[str(n_pep)]
    proteins = proteins[~pd.Series(proteins).isin(bq.index)]

    n_peptides = n_pep
    if n_peptides > top:
        n_peptides = top

    model = bayesquant.compile_model(n_peptides=n_peptides)

#    if bq.shape[0] != 0:
#        human=bq.loc[np.bitwise_and(bq.n_peptides == n_peptides, bq.Organism == "Homo sapiens"),:].shape[0]
#        ecoli=bq.loc[np.bitwise_and(bq.n_peptides == n_peptides, bq.Organism != "Homo sapiens"),:].shape[0]
#    else:
    human=0
    ecoli=0
    print(human)
    print(ecoli)
   
    for p in proteins:
        print(p)
        organism=data.loc[data["protein"] == p, "taxon"].iloc[0]
        if organism=="Homo sapiens":
            human+=1
        else:
            ecoli+=1

        if human>n and organism=="Homo sapiens":
             continue 
        if ecoli>n and organism!="Homo sapiens":
             continue 
        if ecoli>n and human>n:
             break 

        bayesquant.load_data(p, n_peptides)
        trace = bayesquant.fit(model_name=p)
        result = pm.summary(trace, varnames=["estimate"])
        result.index = [p]
        result["Organism"] = organism
        result["n_peptides"] = n_peptides
    
        bq = pd.concat([bq, result])
        bq.to_csv(args["output"], sep="\t")

        bayesquant.traceplot()
        bayesquant.plot_posterior(ref_val=None, xlim=None)
