import pymc3 as pm
from pymc3.variational.callbacks import CheckParametersConvergence
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from tqdm import tqdm
from theano import shared
import shutil
import os
from exp_annot_funcs import *

def protein_model(observed_shared, features_p_shared, x_treat_shared, x_pep_shared, x_run_shared, x_estimate_shared, n_peptides, model_name, n_draws=1000, n_chains=3, hierarchical_center=False, advi=True, remove_backend=True):
 

    """
    Compiles a PYMC3 model to infer the log2FC for protein p
    between two conditions with different treatments,
    i.e the treatment effect

    Keywords:

        -proteins: a string list storing the protein name: ex O43663

        -features: a DataFrame of size nxk where n is the number of peptides available
            in the experiment and k the number of features extracted in a previous step
        
        -sequence: boolean controlling if the sequence features should be used
            to model the variance across peptides of the same protein


        -n_draws: integer. How many draws (sampling) should pymc3 take=?
        -n_chains: integer. How many chains should be run in parallel?

        
   
    The algorithm models the protein quantities with 2 linear regressions:

      The main GLM is used to predict the quant_value from MS1 based on:
          treat: the effect of the treatment
          pep: the effect of the peptide (batch/sequence effect) 
      Other effects may be added, such as fraction
 
      The secondary GLM is built to model pep as a function of features
      extracted from either the peptide or the peptide in combination
      with its surroundings in the protein environment.
      The environment is defined as the 15 residue window on each side
       
      The sequence featues passed should be those listed
      in table 3 from https://www.ncbi.nlm.nih.gov/pubmed/16873510
      Table 3 shows the features that best performed while
      predicting the peptide observability
      
    This function returns a tuple with the trace, the model itself,
    and 2 qc objects
                
   
    """

    if remove_backend and os.path.isdir(model_name):
        shutil.rmtree(model_name)

    # Get the data for the current protein being analyzed
    n_peptides = 2
    n_proteins = 1



    with pm.Model() as model:
            
        intercept = pm.Normal("intercept", 22, 1)
        
        # sigma will be some random error protein specific variance
        sigma = pm.HalfNormal('sigma', 1)
        sigma_pep = pm.HalfNormal('sigma_pep', 1)
        mu_pep = pm.Normal('mu_pep', 0, sigma_pep)
        sigma_treat = pm.HalfNormal('sigma_treat', 1)
        mu_treat = pm.Normal('mu_treat', 0, sigma_treat)
        sigma_run = pm.HalfNormal('sigma_run', 1)
        mu_run = pm.Normal('mu_run', 0, sigma_run)

        if hierarchical_center:
            pep = pm.Normal('pep', mu_pep , sigma_pep, shape = (n_peptides, 1))
            treat = pm.Normal('treat', mu_treat, sigma_treat, shape = (n_proteins*2, 1))
            run = pm.Normal('run', mu_run, sigma_run, shape = (n_proteins*6, 1))
        else:
            pep_offset = pm.Normal("pep_offset", mu=0, sd=1, shape = (n_peptides, 1))
            pep = pm.Deterministic("pep", mu_pep + pep_offset * sigma_pep)
            treat_offset = pm.Normal("treat_offset", mu=0, sd=1, shape=(n_proteins*2, 1))
            treat = pm.Deterministic("treat", mu_treat + treat_offset*sigma_treat)
            run_offset = pm.Normal("run_offset", mu=0, sd=1, shape=(n_proteins*6, 1))
            run = pm.Deterministic("run", mu_run + run_offset*sigma_run)


        # When there is only one protein, thi equals to treat[0] - treat[1]
        estimate = pm.Deterministic('estimate', pm.math.sum(x_estimate_shared.dot(treat), axis=1))
        treatment_effect = pm.Deterministic("treatment_effect", pm.math.sum(x_treat_shared.dot(treat), axis=1))
        peptide_effect = pm.Deterministic("peptide_effect", pm.math.sum(x_pep_shared.dot(pep), axis=1))
        run_effect = pm.Deterministic("run_effect", pm.math.sum(x_run_shared.dot(run), axis=1))

        # Likelihood function for the data.
        mu = pm.Deterministic("mu", intercept + treatment_effect + peptide_effect + run_effect)
        if hierarchical_center:
            obs = pm.Normal("likelihood", mu, sigma, observed=observed_shared)
        else:
            obs_offset = pm.Normal("obs_offset", mu=0, sd=1, shape=(n_peptides*6,1))
            obs = pm.Normal("likelihood", mu+obs_offset*sigma, sigma, observed=observed_shared)


    print("Success: Model {} compiled".format(model_name))

    with model:
        # Parameters of the simulation:
        # Number of iterations and independent chains.
        n_sim = n_draws*n_chains

        # Save traces to the Text backend i.e a folder called model_name containing csv files for each chain
        db = pm.backends.Text('{}'.format(model_name))
        if not advi:
            trace = pm.sample(draws=n_draws, njobs=n_chains, trace=db, tune=2000, nuts_kwargs=dict(target_accept=.95))
        else:
            inference = pm.ADVI()
            approx = pm.fit(n=30000, method=inference, callbacks=[pm.callbacks.CheckParametersConvergence(diff='absolute')])
            trace = approx.sample(draws=n_draws)
            plt.plot(-inference.hist, label='ADVI', alpha=.3)
            plt.legend()
            plt.ylabel('ELBO')
            plt.xlabel('iteration');
            plt.savefig("elbos/ELBO_{}.png".format(model_name))
    
    pm.traceplot(trace, varnames=["estimate"])
    plt.savefig("traceplots/{}.png".format(model_name))
    plt.close()
       
    return model
