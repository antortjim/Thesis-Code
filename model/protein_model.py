import pymc3 as pm
import numpy as np
import matplotlib.pyplot as plt
from theano import shared
import shutil
import os

def protein_model(observed_sh, feats_sh, x_treat_sh, x_pep_sh, x_run_sh, x_estimate_sh,
    n_peptides, model_name, n_draws=1000, n_chains=3,
    hierarchical_center=False, remove_backend=True, sequence=False):

    # Check working environment  
    if not os.path.isdir("traces") or not os.path.isdir("plots/traceplots"):
        msg = "Please create a traces dir and a plots/traceplots dir before running this code"
        raise Exception(msg)

    if remove_backend and os.path.isdir(model_name):
        shutil.rmtree(model_name)

    # The number of proteins in this model is always one
    # i.e this model is fitted protein-wise
    n_prots = 1
    # The number of features is set to 9 for now
    # All peptides have 9 features, stored in feats_shared
    n_features = 9
    

    with pm.Model() as model:
           
        # Build a hyerarchical linear model of
        # log2(MS1 intensities) by accounting for:

             # Peptide effect
             # Run (batch) effect
             # Treatment effect
             # Remaining random effects

        # The difference in treatment effects is an estimate of the log2FC


        # Set a prior on the intercept
        intercept = pm.Normal("intercept", 22, 1)
        
        # Set a prior on the remaining random effects
        sigma = pm.HalfNormal('sigma', 1)

        ## Set priors on the peptide effect
        ################################################
        sigma_pep = pm.HalfNormal('sigma_pep', 1)

        # Not using the sequence
        if not sequence:
            mu_pep = pm.Normal('mu_pep', mu=0, sd=sigma_pep, shape=(n_peptides, 1))

        # Using the peptide sequence
        else: 
            # sequence based modelling
            mu_theta = pm.Normal('theta_generic', 0, sigma_pep, shape = 1)
            theta = pm.Normal('theta', mu_theta, sigma_pep, shape = (n_features, 1))    # 9x1
            theta_inter = pm.Normal('theta_inter', mu_theta, sigma_pep, shape = 1)
            mu_pep = pm.Deterministic("mu_pep", theta_inter + feats_sh.dot(theta)) # n_peptidesx1


        ## Set priors on the treatment and run effects
        ################################################    
        sigma_treat = pm.HalfNormal('sigma_treat', 1)
        mu_treat = pm.Normal('mu_treat', 0, sigma_treat)
        sigma_run = pm.HalfNormal('sigma_run', 1)
        mu_run = pm.Normal('mu_run', 0, sigma_run)

        # Standard implementation of the hyerarchies
        if hierarchical_center:
            pep = pm.Normal("pep", mu_pep, sigma_pep) # n_peptidesx1
            treat = pm.Normal('treat', mu_treat, sigma_treat, shape = (n_prots*2, 1))
            run = pm.Normal('run', mu_run, sigma_run, shape = (n_prots*6, 1))

        # Reparametrization to escape funnel of hell as noted in
        # http://twiecki.github.io/blog/2017/02/08/bayesian-hierchical-non-centered/
        else:
            pep_offset = pm.Normal("pep_offset", mu=0, sd=1, shape = (n_peptides, 1))
            pep = pm.Deterministic("pep", mu_pep + pep_offset * sigma_pep)
            treat_offset = pm.Normal("treat_offset", mu=0, sd=1, shape=(n_prots*2, 1))
            treat = pm.Deterministic("treat", mu_treat + treat_offset*sigma_treat)
            run_offset = pm.Normal("run_offset", mu=0, sd=1, shape=(n_prots*6, 1))
            run = pm.Deterministic("run", mu_run + run_offset*sigma_run)


        # Model the effect for all peptides
        # The sh variables consist of -1,0,1 matrices telling pymc3
        # which parameters shall be used with each peptide
        # In practice, the "clone" each parameter to fit the shape of observed_sh
        # observed_sh is a n_peptides*6x1 tensor
        # The first 6 numbers store the MS1 intensities of the first peptide in the 6 runs
        # The next 6 those of the second peptide, and so on

        estimate = pm.Deterministic('estimate', pm.math.sum(x_estimate_sh.dot(treat), axis=1))
        treatment_effect = pm.Deterministic("treatment_effect", pm.math.sum(x_treat_sh.dot(treat), axis=1))
        peptide_effect = pm.Deterministic("peptide_effect", pm.math.sum(x_pep_sh.dot(pep), axis=1))
        run_effect = pm.Deterministic("run_effect", pm.math.sum(x_run_sh.dot(run), axis=1))

        # BIND MODEL TO DATA
        mu = pm.Deterministic("mu", 
            intercept + treatment_effect + peptide_effect + run_effect) #n_peptides*6x1
        if hierarchical_center:
            obs = pm.Normal("obs", mu, sigma, observed=observed_sh)
        else:
            obs_offset = pm.Normal("obs_offset", mu=0, sd=1, shape=(n_peptides*6,1))
            obs = pm.Normal("obs", mu+obs_offset*sigma, sigma, observed=observed_sh)


    print("Success: Model {} compiled".format(model_name))

    with model:
        # Parameters of the simulation:
        # Number of iterations and independent chains.
        n_sim = n_draws*n_chains

        # Save traces to the Text backend i.e a folder called
        # model_name containing csv files for each chain
        trace_name = 'traces/{}'.format(model_name)
        db = pm.backends.Text(trace_name)
        trace = pm.sample(draws=n_draws, njobs=n_chains, trace=db,
                          tune=2000, nuts_kwargs=dict(target_accept=.95))
    
    # Save a traceplot 
    pm.traceplot(trace, varnames=["estimate"])
    traceplot = "plots/traceplots/{}.png".format(model_name)
    plt.savefig(traceplot)
    plt.close()
       
    return model