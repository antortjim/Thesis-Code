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
from MSBay import MSBay
n=20
top=20

data = pd.read_csv("data/ms1_intensities.tsv", sep = "\t")
data['H1'] = data['H1'].astype(theano.config.floatX)
data['H2'] = data['H2'].astype(theano.config.floatX)
data['H3'] = data['H3'].astype(theano.config.floatX)
data['L1'] = data['L1'].astype(theano.config.floatX)
data['L2'] = data['L2'].astype(theano.config.floatX)
data['L3'] = data['L3'].astype(theano.config.floatX)


counts0 = np.unique(data.protein, return_counts=True)
protein_counts = {str(c): counts0[0][counts0[1] == c] for c in np.unique(counts0[1])}
MSBayQ_file="data/MSBayQ.tsv"
if os.path.isfile(MSBayQ_file):
    MSBayQ = pd.read_csv(MSBayQ_file, sep="\t", index_col=0)
else:
    MSBayQ = pd.DataFrame({"mean": [], "sd": [], "hpd_2.5":[], "hpd_97.5":[],"Organism":[], "n_peptides":[]})


for n_pep in range(2,6):
    proteins = protein_counts[str(n_pep)]
    proteins = proteins[~pd.Series(proteins).isin(MSBayQ.index)]
    msbay = MSBay(data.loc[data["protein"].isin(proteins)], features=None)

    n_peptides = n_pep
    if n_peptides > top:
        n_peptides = top

    model = msbay.compile_model(n_peptides=n_peptides)

    if MSBayQ.shape[0] != 0:
        human=MSBayQ.loc[np.bitwise_and(MSBayQ.n_peptides == n_peptides, MSBayQ.Organism == "Homo sapiens"),:].shape[0]
        ecoli=MSBayQ.loc[np.bitwise_and(MSBayQ.n_peptides == n_peptides, MSBayQ.Organism != "Homo sapiens"),:].shape[0]
    else:
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

        msbay.load_data(p, n_peptides)
        trace = msbay.fit(model_name=p)
        result = pm.summary(trace, varnames=["estimate"])
        result.index = [p]
        result["Organism"] = organism
        result["n_peptides"] = n_peptides
    
        MSBayQ = pd.concat([MSBayQ, result])
        MSBayQ.to_csv("data/MSBayQ.tsv", sep="\t")

