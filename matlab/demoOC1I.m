% DEMOC1I Demo of probabilistic PCA with noise on OC1I dataset. The dataset contains ~13K genes and 12 samples; approximate running time 18hrs.

% NPPCA

% Fix a seed so that results are repeatable.
randn('seed', 2e5);
rand('seed', 2e5);
colordef white
% We request four latent dimensions, but model only uses 3.
latentDim = 7;

% Set up the options.
options = nppcaOptions;

% Load first 20 points from OC1 data.
[probes, annotations, Y, varY] = nppcaLoadData('OC1');

[model, expectations] = nppcaMaster(Y, varY, latentDim, options);
