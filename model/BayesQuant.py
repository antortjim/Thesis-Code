import pymc3 as pm
import numpy as np
import matplotlib.pyplot as plt
from theano import shared
import theano
import shutil
import os
from exp_annot_funcs import *
import sys
import pickle
from sklearn.decomposition import PCA 
import seaborn as sns


class BayesQuant:

    "Compile a BayesQuant model, load data into it and posterior distributions"

    def __init__(self, data=None, features=None, plot_dir="../../Report/plots"):
        self.data = data
        self.features = features
        self.p = None
        self.trace = None
        self.pca = None
        self.plot_dir = plot_dir


    def read_data(self, data_path="data/data.tsv", features_path="data/advanced_features.tsv"):
        data = pd.read_csv(data_path, sep = "\t")
        for c in data.columns[2:].values:
            data[c] = data[c].astype(theano.config.floatX)


        self.data = data


        if features_path is not None:
            features = pd.read_csv(features_path, sep = "\t").iloc[:,1:]
            features=features.iloc[data.index,:]
            n_seq_feat = features.shape[1]
            print("Number of sequence features {}".format(n_seq_feat))
            n_prots = len(np.unique(data.protein))
            print("Number of proteins {}".format(n_prots))
            n_peptides = data.shape[0]
            print("Number of peptides {}".format(n_peptides))
            print("This is the first row of the features array:")
            print(features.shape)
            print(features.iloc[:1,:])

            self.features = features

        return (self.data, self.features)




    def compile_model(self, n_peptides, hierarchical_center=False):


        if self.features is None:
             sequence=False
        elif isinstance(self.features, pd.DataFrame):
             sequence = True
        else:
             sys.exit("Please set features to either None or a pd.DataFrame")


        self.observed_sh  = shared(np.array([0.,]*6*n_peptides))
        self.feats_sh  = shared(np.array([[0.,]*9,]*n_peptides))
        self.x_treat_sh  = shared(np.array([[0.,]*2,]*6*n_peptides))
        self.x_pep_sh  = shared(np.array([[0.,]*n_peptides,]*6*n_peptides))
        self.x_run_sh  = shared(np.array([[0.,]*6,]*6*n_peptides))
        self.x_estimate_sh  = shared(np.array([[0.,]*2,]*1))
     
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
    
    
            # Set a prior on the random effects
            #sigma = pm.HalfNormal('sigma', 1)
            # Set a prior on the intercept
            intercept = pm.Normal("intercept", mu=22, sd=1)
            
    
            ## Set priors on the peptide effect
            ################################################
            sigma_pep = pm.HalfNormal('sigma_pep', 1)
    
            # Not using the sequence
            if not sequence:
                mu_pep = pm.Normal('mu_pep', mu=0, sd=sigma_pep)
    
            # Using the peptide sequence
            else: 
                # sequence based modelling
                mu_theta = pm.Normal('theta_generic', 0, sigma_pep)
                theta = pm.Normal('theta', mu_theta, sigma_pep, shape = (n_features, 1))    # 9x1
                theta_inter = pm.Normal('theta_inter', mu_theta, sigma_pep, shape = 1)
                mu_pep = pm.Deterministic("mu_pep", theta_inter + self.feats_sh.dot(theta))
    
    
            ## Set priors on the treatment and run effects
            ################################################    
            sigma_treat = pm.HalfNormal('sigma_treat', 1)
            mu_treat = pm.Normal('mu_treat', 0, sigma_treat)
            sigma_run = pm.HalfNormal('sigma_run', 1)
            mu_run = pm.Normal('mu_run', 0, sigma_run)
    
            # Standard implementation of the hyerarchies
            if hierarchical_center:
                pep = pm.Normal("pep", mu_pep, sigma_pep, shape = (n_peptides, 1)) # n_peptidesx1
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
    
            #estimate = pm.Deterministic('estimate', pm.math.sum(self.x_estimate_sh.dot(treat), axis=1))
            estimate = pm.Deterministic('estimate', self.x_estimate_sh.dot(treat))
            treatment_effect = pm.Deterministic("treatment_effect", pm.math.sum(self.x_treat_sh.dot(treat), axis=1))
            peptide_effect = pm.Deterministic("peptide_effect", pm.math.sum(self.x_pep_sh.dot(pep), axis=1))
            run_effect = pm.Deterministic("run_effect", pm.math.sum(self.x_run_sh.dot(run), axis=1))
    
            # BIND MODEL TO DATA
            mu = intercept + treatment_effect + peptide_effect + run_effect
            #mu = pm.Deterministic("mu", effects) #n_peptides*6x1
            epsilon = pm.HalfNormal('epsilon', sd=1)
            if hierarchical_center:
                obs = pm.Normal("obs", mu, epsilon, observed=self.observed_sh)
            else:
                obs_offset = pm.Normal("obs_offset", mu=0, sd=1, shape=(n_peptides*6,1))
                obs = pm.Normal("obs", mu+obs_offset*epsilon, epsilon, observed=self.observed_sh)
    
    
        print("Success: Model compiled")

        self.model = model
        return model
    
    def sample(self, model_name=None, n_draws=1000, n_chains=3, remove_backend=True):

        p = self.p
        if model_name is not None:
            p = model_name    

        # Check working environment  
        if not os.path.isdir("traces") or not os.path.isdir("plots/traceplots"):
            msg = "Please create a traces dir and a plots/traceplots dir before running this code"
            raise Exception(msg)
    
        if remove_backend and os.path.isdir("traces/{}".format(p)):
            shutil.rmtree("traces/{}".format(p))
    
    
        with self.model:
            # Parameters of the simulation:
            # Number of iterations and independent chains.
            n_sim = n_draws*n_chains
    
            # Save traces to the Text backend i.e a folder called
            # model_name containing csv files for each chain
            trace_name = 'traces/{}'.format(p)
            db = pm.backends.Text(trace_name)
            trace = pm.sample(draws=n_draws, njobs=n_chains, trace=db,
                              tune=2000, nuts_kwargs=dict(target_accept=.95))
        
        # Save a traceplot 
        pm.traceplot(trace, varnames=["estimate"])
        traceplot = "plots/traceplots/{}.png".format(p)
        plt.savefig(traceplot)
        plt.close()


        self.trace = trace
           
        return trace

    def fit(self, model_name=None, n_draws=40000):

        p = self.p
        if model_name is not None:
            p = model_name
 
        try:
            os.mkdir("traces/{}".format(p))
        except:
            print("Dir exists")


        with self.model:
            
            inference = pm.ADVI()
            # how can the trace be saved when using pm.fit??
            trace = inference.fit(n=n_draws).sample()



            with open("traces/{}/trace.pik".format(p), 'wb') as f:
                pickle.dump({'model': self.model, 'trace': trace}, f)

            #with open('trace.p', 'rb') as f:
            #    test1 = pickle.load(f)

            # trace = pm.fit(n=n_draws, method=inference).sample()


        plt.plot(-inference.hist, alpha=.5)
        plt.legend()
        plt.ylabel('ELBO')
        plt.xlabel('iteration');
        plt.savefig("plots/ELBO/{}".format(p))
        plt.close()


 
        df=pd.DataFrame({"estimate": trace["estimate"][:,0,0]})
        df.to_csv("traces/{}/chain-0.tsv".format(p))

        self.trace = trace
        
        return trace
    
    def load_data(self, p, top=3):
    
        df = pd.DataFrame({"std": np.std(self.data.loc[self.data.protein == p,:].iloc[:,2:5].values, axis=1) + np.std(self.data.loc[self.data.protein == p,:].iloc[:,5:8].values, axis=1)})
        best_data = self.data.loc[self.data.protein == p,:].iloc[df.sort_values(by="std", ascending=True).iloc[:top,:].index,:]
        variables = create_variables(best_data, self.features, [p])
        observed, feats, x_treat, x_pep, x_run, x_estimate = variables

        self.observed_sh.set_value(observed)
        if feats is not None:
            self.feats_sh.set_value(feats)    

        self.x_treat_sh.set_value(x_treat)
        self.x_pep_sh.set_value(x_pep)
        self.x_run_sh.set_value(x_run)
        self.x_estimate_sh.set_value(x_estimate)
        self.p = p


    def ppc(self, samples=500):

        with self.model:
            sim = pm.sample_ppc(self.trace, samples=samples)["obs"]

        sim=sim[:,0,:]
        peps = list(map(lambda i: sim[:, (i*6):((i*6)+6)], range(n_peptides)))
        estimates = list(map(lambda x: np.mean(x[:,:3], axis=1) - np.mean(x[:,3:6], axis=1), peps))
        fig, ax = plt.subplots(1,n_peptides, figsize=(15,5))
        [ax[i].hist(e) for i, e in enumerate(estimates)]
        plt.savefig(os.path.join(self.plot_dir, "PPC/histogram_{}".format(p)))
        plt.close()


        peps_data = np.vstack(peps)
        sim_transformed = self.pca.transform(peps_data)
        sim_transformed
        plt.scatter(sim_transformed[:,0], sim_transformed[:,1], label="Sim obs")
        obs = self.data.loc[self.data.protein == p].iloc[:,2:]
        obs_transformed = pca.transform(obs)
        plt.scatter(obs_transformed[:,0], obs_transformed[:,1], c = "red", label="True obs")
        plt.legend()
        plt.xlabel("PC1")
        plt.ylabel("PC2")
        plt.savefig(os.path.join(self.plot_dir, "PPC/PCA_{}".format(p)))
        plt.close()

    def PCA(self):
        x = self.data.values[:,2:]
        pca = PCA(n_components=2)
        pca.fit(x)
        self.pca = pca
        print("Percentage of variance explained by PC1 and PC2")
        print(pca.explained_variance_ratio_*100)
        x_transformed = pca.transform(x)
        pca_data = pd.DataFrame({"taxon": self.data.taxon, "PC1": x_transformed[:,0], "PC2": x_transformed[:,1]})
        myPlot = sns.FacetGrid(col="taxon", hue='taxon', data=pca_data, size=5)
        myPlot = myPlot.map(plt.scatter, "PC1", "PC2", alpha=0.3)
        myPlot = myPlot.map_dataframe(plt.plot, [min(pca_data.PC1),max(pca_data.PC1)], [0, 0], 'r-').add_legend().set_axis_labels("PC1", "PC2")
        plt.savefig(os.path.join(self.plot_dir, "PCA.png"))