MS-Bay - Bayesian estimation of log2(FC)
============================================


 * [Quickstart](#quickstart)
 * [Introduction](#introduction)
 * [Quantification](#quantification)
 * [Implementation](#implementation)

---


## Quickstart



###  Read data for a full protemics experiment dataset

```
# data.tsv was created running generate_model_dataset.R using the peptides.txt and proteinGroups.txt files from the MaxLFQ proteome benchmark dataset
# It was processed using MSqRob::preprocess_MaxQuant()
data = pd.read_csv("data/data.tsv", sep = "\t")
data.head(10)
```
![skema](figures/MS1_intensities.png)


#### Instantiate an MsBay object
```
from MSBay import MSBay
msbay = MSBay(data, features)
```

#### Compile a model for 2 (n) peptides
```
model = msbay.compile_model(n_peptides=2)
```

#### Load the peptides for protein P76397.

Its true fold change is log2(3) = 1.58
```
# An E. coli protein with 2 peptides observed.
msbay.load_data("P76397")
```

#### (A) Compute posterior by using the NUTS sampler

The NUTS sampler takes longer but it's guaranteed to find the true posterior at convergence, i.e given enough samples.
```
trace_nuts = msbay.sample(model_name=p,n_draws=1000, n_chains=3)
```
```
Auto-assigning NUTS sampler...
Initializing NUTS using jitter+adapt_diag...
/home/antortjim/anaconda3/envs/bayesian/lib/python3.6/site-packages/pymc3/model.py:384: FutureWarning: Conversion of the second argument of issubdtype from `float` to `np.floating` is deprecated. In future, it will be treated as `np.float64 == np.dtype(float).type`.
  if not np.issubdtype(var.dtype, float):
Multiprocess sampling (3 chains in 3 jobs)
NUTS: [obs_offset, run_offset, treat_offset, pep_offset, mu_run, sigma_run_log__, mu_treat, sigma_treat_log__, mu_pep, sigma_pep_log__, sigma_log__, intercept]
100%|██████████| 3000/3000 [01:56<00:00, 25.77it/s]
There were 24 divergences after tuning. Increase `target_accept` or reparameterize.
There were 11 divergences after tuning. Increase `target_accept` or reparameterize.
There were 17 divergences after tuning. Increase `target_accept` or reparameterize.
The number of effective samples is smaller than 25% for some parameters.
```


#### (B) Compute posterior by using an approximation method (ADVI)

The ADVI engine takes very little and finds a good approximation to the posterior.
```
trace_advi = msbay.fit(model_name=p, n_draws=40000)
```

```
Average Loss = 192.22: 100%|██████████| 40000/40000 [00:27<00:00, 1453.28it/s]
Finished [100%]: Average Loss = 192.24
No handles with labels found to put in legend.
```

#### Plot the posterior probability distribution for the log2FC

**NUTS**
```
pm.plot_posterior(trace_nuts, varnames=["estimate"], ref_val=np.log2(3), color='LightSeaGreen', rope=[-0.4, 0.4])
plt.savefig("plots/posteriors/{}_nuts.png".format(p))
```

**ADVI**
```
pm.plot_posterior(trace_advi, varnames=["estimate"], ref_val=np.log2(3), color='LightSeaGreen', rope=[-0.4, 0.4])
plt.savefig("plots/posteriors/{}_advi.png".format(p))
``` 
![](plots/posteriors/P76397_nuts.png)
![](plots/posteriors/P76397_advi.png)


Both methods return a posterior distribution whose 95% HDI contain the true parameter (log2(3) = 1.58). However, the ADVI method, running Variational Inference, proofs to be much faster and much more accurate than the standard NUTS sampler. The presence of several layers or hyerarchies in the model gives NUTS a hard time. ADVI gets over these problems by performing an approximation that nevertheless returns very good results.

As a consequence, the NUTS HDI partially overlaps the ROPE (region of practical equivalence), while the ADVI HDI clearly does not.

If we use the ADVI engine to quantify 10 new proteins, we get the following results:

![](plots/performance.png)

The posterior  95% HDI overlaps a ROPE defined as -.4,.4 for all tested *Homo sapiens* proteins, while it does not for *E. coli* proteins, which cluster around the log2(3), the true log2FC for them.

## Introduction  

Protein quantification in proteomics is used to better discern how a biological system responds to a specific stimuli at the protein level. It enables the extraction of quantitative information on its protein composition, which could lead to a deeper understanding of many cell processes in which proteins are involved.


Quantification can be achieved by using MS1 intensity measurements as proxy for protein quantities. However, in order to correctly asses the treatment effect i.e the changes in protein quantities due to the actual treatments applied, one needs to factor in the peptide variability and the run (batch) variability.

This module performs label-free relative protein quantification using a simple Bayesian model while accounting for these effects:

  * treat: the effect of the treatment.
  * pep: the effect of the peptide (optionally modeled from sequence).
  * run: the effect of the run (batch).


It implements a linear regression model similar to that from the [MSqRob package](https://github.com/statOmics/MSqRob). While the original MSqRob minimises an OLS loss function with ridge regression and Huber weights, the present program instead runs Monte Carlo-Markov Chain (MCMC) simulations to determine the posterior probability distribution of the model parameters. Once the model parameters are fitted to the data for one protein, the posterior probability of its fold change estimate can be infered.

## Quantification

The tool is a new alternative to the "Protein quantification" step. Current alternatives are 

* MSnBase: An R Bioconductor package performing summarization methods
* MSqRob: An R package implementing a peptide-based model. Supports moFF, MaxQuant and open generic formats.
* MaxLFQ: Part of the MaxQuant suite, performs peptide summarization and least squares regression.
* StPeter: Part of the TPP suite, implements the spectral index method. Only TPP input formats are supported.
* Non-free software: Progenesis QI, Proteome Discoverer, GeneData, etc.

Like MSqRob, this package performs relative quantification, i.e it provides an estimate of the abundance ratio between 2 conditions, but it does not provide an absolute measurement of the quantities on each one. Thus, its results can be regarded as an estimate of the log2FC, or the log2 abundance ratio between the two conditions being compared.

![skema](figures/proteomics_skema.png)


## Implementation


GRAPH_OF_MODEL



The code for the model is contained in `protein_model.py`. It contains the `compute_posterior()` method:

```
def compute_posterior(self, n_peptides, hierarchical_center=False, remove_backend=True, sequence=False):

    """
    Keywords
    
    -shared_variables: np.arrays made into Tensors with theano.shared() storing
        all the data for a single protein
        The data contained in them can be changed with var.set_value()
        for fast fitting of the model to several protens.

    -n_peptides: int, number of peptides modelled. Changes in this
        argument define different models that can be fitted
        to proteins with different number of peptides

    -model_name: string coding the name of the trace and the traceplot

    -n_draws: int. How many draws (sampling) should pymc3 take=?
    -n_chains: int. How many chains should be run in parallel?

    -hierarchical_center: if False, perform reparametrization trick from
    
    -remove_backend=if True, overwrite an existing backend if it exists
    
    -sequence: boolean controlling if the sequence features should be used
        to model the variance across peptides of the same protein
    """
```

The algorithm models the log2fc estimate with 2 linear regressions:

The main GLM is used to predict the quant_value from MS1 based on the effects mentioned above. Other effects may be added, such as fraction.

The secondary GLM, if sequence=True, is built to model pep as a
function of features extracted from either the peptide sequence alone or in combination
with its surroundings in the protein environment.
The environment is defined as the 15 residue window on each side.
The result of this LM is an estimate of the peptide bias for each peptide (2,3,4..)
based on the features

If sequences=False, this estimate is modelled with a simple Normal distribution.

The sequence features passed should be those listed in table 3 from https://www.ncbi.nlm.nih.gov/pubmed/16873510
Table 3 shows the features that best performed while predicting the peptide observability.

This function returns the model object.