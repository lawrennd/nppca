% DEMOC1p27 Demo of probabilistic PCA with noise on reduced OC1 dataset. We consider only 500 genes including the p27
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
probes=probes(3501:4000);
annotation=annotation(3501:4000);
Y=Y(3501:4000,:);
varY=varY(3501:4000,:)/10;

[model, expectations] = nppcaMaster(Y, varY, latentDim, options); 

save risultatip27 model expectations