import numpy as np
import pandas as pd
from tqdm import tqdm

def create_x_pep(data):

    observed=data.values[:,2:].flatten().astype(float)
    x_pep = np.zeros((len(observed), data.shape[0]))
    k=0
    for i in range(x_pep.shape[1]):
        x_pep[(0+k):(6+k),i] = 1
        k+=6


    print("Shape of x_pep is {}, {}".format(*x_pep.shape))
    return x_pep

def create_x_treat(data, n_proteins):
    observed=data.values[:,2:].flatten().astype(float)
    x_treat = np.zeros((len(observed), n_proteins * 2))
    
    proteins_column = data["protein"].values.tolist()
    last_protein = proteins_column[0]
    colIndex=0
    rowIndex=0
    acum = 0
    for i, p in tqdm(enumerate(proteins_column)):
        if p != last_protein:
            colIndex += 2
            last_protein=p
            acum += 1
        x_treat[(0+rowIndex):(3+rowIndex),colIndex] = 1
        x_treat[(3+rowIndex):(6+rowIndex),colIndex+1] = 1
    
        rowIndex += 6

    print("Shape of x_treat is {}, {}".format(*x_treat.shape))
    return x_treat

def create_x_run(data):
    observed=data.values[:,2:].flatten().astype(float)
    n_proteins = len(np.unique(data.protein))
    n_runs = 6
    x_run = np.zeros((len(observed),n_proteins*n_runs))

    last_protein = data.protein.values.tolist()[0]
    colIndex = 0
    rowIndex = 0
    acum = 0
    for p in data["protein"].values.tolist():
        if p != last_protein:
            colIndex += n_runs
            last_protein=p
            acum += 1

        for offset in range(n_runs):
            x_run[(rowIndex+offset), (colIndex+offset)]=1
        
        rowIndex += n_runs

    return x_run

         
def create_x_estimate(n_proteins):
    x_estimate = np.zeros((n_proteins, n_proteins*2))
    i=0
    j = 0
    for i in range(x_estimate.shape[0]):
        x_estimate[i,j]=1
        j+=1 
        x_estimate[i,j]=-1
        j+=1

    print("Shape of x_estimate is {}, {}".format(*x_estimate.shape))
    return x_estimate
    

def create_variables(data_p, features, proteins):
    indices = data_p.index
    if features is None:
        features_p = None
    else:
        features_p = features.iloc[indices,:].values
        n_features = features_p.shape[1]

    n_peptides = len(indices)
    n_proteins = len(proteins)

    observed=data_p.values[:,2:].flatten().astype(float)
    #log2fc = np.array([0,0,0, log2fc, log2fc, log2fc,]*n_peptides)
    x_treat = create_x_treat(data_p, n_proteins)
    x_pep = create_x_pep(data_p)
    x_run = create_x_run(data_p)
    x_estimate = create_x_estimate(n_proteins)

    return (observed, features_p, x_treat, x_pep, x_run, x_estimate)

def load_stats(data_p, features, proteins):
    data_p=data.loc[data.protein.isin(proteins),:]
    
    observed=data_p.values[:,2:]
    mean_intercept = np.mean(observed[:,:3].flatten())
    observed=observed.flatten().astype(float)
    observed = np.expand_dims(observed, axis=0)

    indices = data_p.index
    features_p = features.iloc[indices,:].values
    n_peptides = len(indices)
    n_proteins = len(proteins)
    n_features = features_p.shape[1]


    return (n_peptides, n_proteins, n_features, mean_intercept)

