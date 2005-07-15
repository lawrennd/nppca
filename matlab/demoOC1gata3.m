% DEMOOC1GATA3 Simple demo of probabilistic PCA with noise on reduced OC1 dataset. 

% NPPCA

% We consider only 500 genes including the gata3
% transcription factor, whose profile is then reconstructed using the
% denoising process.

% Fix a seed so that results are repeatable.
randn('seed', 2e5);
rand('seed', 2e5);
colordef white
% We request four latent dimensions, but model only uses 3.
latentDim = 7;

% Set up the options.
options = nppcaOptions;

% Load first 20 points from OC1 data.
[probes,annotation,Y,varY] = nppcaLoadData('OC1B');
probes=probes(5501:6000);
annotation=annotation(5501:6000);
Y=Y(5501:6000,:);
varY=varY(5501:6000,:);

[model, expectations] = nppcaMaster(Y, varY, latentDim, options);
save resultsGata3 model expectations Y varY
figure, nppcaProfilePlotter(model,expectations,Y, varY,454,'gata3')