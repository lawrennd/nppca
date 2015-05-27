NPPCA software
Version 0.111		Saturday 25 Feb 2006 at 10:16
Copyright (c) 2006 Guido Sanguinetti and Neil D. Lawrence

This is the software associated with the Bioinformatics paper:

G. Sanguinetti, M. Milo, M. Rattray and N. D. Lawrence.  (2005) "Accounting for probe-level noise in principal component analysis of microarray data" in Bionformatics 21 (19), pp 3748--3754

The software implements an E-M algorithm for recovering the principal components in the presence of noise. 

Version 0.111
-------------

Update containing a missing file nppcaLikelihoodCheck.m

MATLAB Files
------------

Matlab files associated with the toolbox are:

clusterLoadData.m: Load a dataset for demonstrating noisy PPCA.
demCNS.m: Simple demo of probabilistic PCA with noise on brain tumour dataset.
demNppcaLoadData.m: Load a dataset for demonstrating noisy PPCA.
demoOC1gata3.m: Simple demo of probabilistic PCA with noise on reduced OC1 dataset. 
demoOC1I.m: Demo of probabilistic PCA with noise on OC1I dataset. The dataset contains ~13K genes and 12 samples; approximate running time 18hrs.
gmosReadTxt.m: reads TXT file for the OC1 data files.
nppcaEstep.m: Estep for the noisy probabilistic PCA.
nppcaFAInit.m: Initialises NPPCA with factor analysis with constant errors
nppcaInit.m: Initialise a Noisy probabilistic PCA model.
nppcaInit100.m: Initialise a Noisy probabilistic PCA model.
nppcaLikelihoodBound.m: Gives the variational bound on the log likelihood for noisy PPCA.
nppcaLikelihoodCheck.m: Compute the difference in likelhoods.
nppcaLoadData.m: Load a dataset for demonstrating noisy PPCA.
nppcaMaster.m: performs the noisy probabilistic PCA algorithm on a dataset.
nppcaNewtonUpdateLogSigma.m: performs a Newton update for the logsigma dependence of the likelihood.
nppcaNewtonUpdateSigma.m: performs a Newton update for the sigma dependence of the likelihood.
nppcaOptions.m: Set up an options vector for noisy PPCA.
nppcaProfilePlotter.m: shows comparison of reconstructed and original profile for genes.
nppcaRemoveRedundancy.m: Remove the redundancy in m and Cinv.
nppcaSigmaGradient.m: Wrapper function for gradient with respect to sigma.
nppcaSigmaObjective.m: Wrapper function for objective as a function of Sigma.
nppcaUpdateCinv.m: Update the latent precision for the noisy PPCA model.
nppcaUpdateM.m: Update the latent mean for the noisy PPCA model.
nppcaUpdateMu.m: Update the mean for the noisy PPCA model.
nppcaUpdateW.m: Update the W matrice for a noisy probabilistic PCA model.
RALoadData.m: Load a dataset for demonstrating noisy PPCA.
reconstruct.m: reproduces a denoised reconstruction of the gene expressions with errorbars.
