Bayesian estimation of log2FC
=========================================0

This module performs label-free relative protein quantification using a simple Bayesian model


Compiles a PYMC3 model to infer the log2FC for protein p
between two conditions with different treatments,
i.e the treatment effect

Keywords:

    -shared_variables: np.arrays made into Tensors with theano.shared() storing
        all the data for a single protein
        The data contained in them can be changed with var.set_value()
        for fast fitting of the model to several protens.

    -n_peptides: number of peptides modelled. Changes in this
        argument define different models that can be fitted
        to proteins with different number of peptides

    -model_name: string coding the name of the trace and the traceplot

    -n_draws: integer. How many draws (sampling) should pymc3 take=?
    -n_chains: integer. How many chains should be run in parallel?

    -hierarchical_center: if False, perform reparametrization trick from
    
    -remove_backend=if True, overwrite the an existing backend if it exists
    
    -sequence: boolean controlling if the sequence features should be used
        to model the variance across peptides of the same protein      

The algorithm models the log2fc estimate with 2 linear regressions:

  The main GLM is used to predict the quant_value from MS1 based on:
      treat: the effect of the treatment
      pep: the effect of the peptide (sequence effect)
      run: the effect of the run (batch)

  Other effects may be added, such as fraction

  The secondary GLM, if sequence=True, is built to model pep as a
  function of features extracted from either the peptide or the peptide in combination
  with its surroundings in the protein environment.
  The environment is defined as the 15 residue window on each side.
  The result of this LM is an estimate of the peptide bias for each peptide (2,3,4..)
  based on the features

  If sequences=False, this estimate is modelled with a simple Normal distribution.
   
  The sequence features passed should be those listed
  in table 3 from https://www.ncbi.nlm.nih.gov/pubmed/16873510
  Table 3 shows the features that best performed while
  predicting the peptide observability.

  
This function returns the model object
