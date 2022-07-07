# Bayesian-covariance-estimation
 
Simulation experiments code for comparing the following approaches for covariance matrix estimation are provided. 

SCOV: Basic sample covariance

GLasso: A Frequentist approach. Sparse inverse covariance estimation, Friedman et al 2007. https://tibshirani.su.domains/ftp/graph.pdf

SBIFM: Sparse Bayesian infinite factor models, Bhattacharya & Dunson 2011, https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3419391/pdf/asr013.pdf

STEP: Precision matrix estimation. Simultaneous Regression and Covariance Selection, Samanta et al. 2022. https://arxiv.org/abs/2201.05653

BGLasso: Precision matrix estimation. Bayesian Graphical Lasso, Wang 2012. https://projecteuclid.org/journals/bayesian-analysis/volume-7/issue-4/Bayesian-Graphical-Lasso-Models-and-Efficient-Posterior-Computation/10.1214/12-BA729.full


----------------------------------------------

code.R contains the function to run during the simulation

utils.R contains useful functions for calculating the performance measures

simulation.R contains code for the data generation, and running the simulation experiments. 
