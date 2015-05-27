ReadMe file for the NPPCA toolbox version 0.11 Monday, Jul 18, 2005 at 14:25:02
Written by Guido Sanguinetti and Neil D. Lawrence.

License Info
------------

This software is free for academic use. Please contact Neil Lawrence if you are interested in using the software for commercial purposes.

This software must not be distributed or modified without prior permission of the author.



File Listing
------------

RALoadData.m: Load a dataset for demonstrating noisy PPCA.
clusterLoadData.m: Load a dataset for demonstrating noisy PPCA.
demCNS.m: Simple demo of probabilistic PCA with noise on brain tumour dataset.
demNppcaLoadData.m: Load a dataset for demonstrating noisy PPCA.
demoOC1I.m: Demo of probabilistic PCA with noise on OC1I dataset. The dataset contains ~13K genes and 12 samples; approximate running time 18hrs.
demoOC1gata3.m: Simple demo of probabilistic PCA with noise on reduced OC1 dataset. 
gmosReadTxt.m: reads TXT file for the OC1 data files.
nppcaEstep.m: Estep for the noisy probabilistic PCA.
nppcaFAInit.m: Initialises NPPCA with factor analysis with constant errors
nppcaInit.m: Initialise a Noisy probabilistic PCA model.
nppcaInit100.m: Initialise a Noisy probabilistic PCA model.
nppcaLikelihoodBound.m: Gives the variational bound on the log likelihood for noisy PPCA.
nppcaLoadData.m: Load a dataset for demonstrating noisy PPCA.
nppcaMaster.m: performs the noisy probabilistic PCA algorithm on a dataset.
nppcaNewtonUpdateLogSigma.m: performs a Newton update for the logsigma dependence of the likelihood.
nppcaNewtonUpdateSigma.m: performs a Newton update for the sigma dependence of the likelihood.
nppcaOptions.m: Set up an options vector for noisy PPCA.
nppcaRemoveRedundancy.m: Remove the redundancy in m and Cinv.
nppcaSigmaObjective.m: Wrapper function for bjective as a function of Sigma.
nppcaUpdateCinv.m: Update the latent precision for the noisy PPCA model.
nppcaUpdateM.m: Update the latent mean for the noisy PPCA model.
nppcaUpdateMu.m: Update the mean for the noisy PPCA model.
nppcaUpdateW.m: Update the W matrice for a noisy probabilistic PCA model.
reconstruct.m: reproduces a denoised reconstruction of the gene expressions with errorbars.
