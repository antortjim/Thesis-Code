# Thesis-Code

This repository stores the code written by Antonio Ortega during his Master Thesis in Bioinformatics at the University of Copenhagen in Denmark (Københavns Universitet). The code is split in 2 folders:

* [`pipeline`](https://github.com/antortjim/Thesis-Code/tree/master/pipeline): bash, R and Python scripts required to run an OS label-free proteomics quantification pipeline in Linux using the tools published by the [Compomics group at University of Ghent](https://compomics.com/).

* [`model`](https://github.com/antortjim/Thesis-Code/tree/master/model): Python scripts powering BayesQuant, a program that uses probabilistic programming to compute relative quantities from MS1 intensity data, in a format similar to that used by [MSqRob](https://github.com/statOmics/MSqRob). It consists of a PyMC3 model.

Additionally, the thp1 folder stores data in different processing stages analysed by the pipeline for a benchmark dataset of THP1 cells.

[Another repo](https://github.com/antortjim/Thesis-Report)  contains the latex code required to compile the PDF document, available [here](http://people.binf.ku.dk/rnq313/master_thesis/thesis.pdf) .
  