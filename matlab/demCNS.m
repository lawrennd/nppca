% DEMCNS Simple demo of probabilistic PCA with noise on brain tumour dataset.

% NPPCA

% Fix a seed so that results are repeatable.
randn('seed', 2e5);
rand('seed', 2e5);

% Request 12 latent dimensions.
latentDim = 12;

% Set up the options.
options = nppcaOptions;

% Load first 20 points from OC1 data.
[probeNames, Y] = readCNS_DATA('../data/brainData_exprs.csv');
[probeNames, varY] = readCNS_DATA('../data/brainData_standDev.txt');

varY=varY.^2;
% set zeros to 1e-6;
varY(find(varY<1e-6))=1e-6;
dataDim = size(Y, 2);
numData = size(Y, 1);

% Initialise the model --- reset to PCA.
[model, expectations] = nppcaInit(Y, varY, latentDim);

% Compute the starting likelihood.
maxDeltaL = 1;
counter=0;
model.mu = zeros(size(model.mu));
params = {'mu', 'W', 'sigma', 'estep'};
oldL = nppcaLikelihoodBound(model, expectations, varY, Y);

while (maxDeltaL > options.tol & counter < options.maxIters)
  maxDeltaL = 0;
  counter=counter+1;
  if counter == 1
    order = randperm(length(params));
  end
  for i = order
    updateParam = params{i};
    switch updateParam
     case 'mu'
      model.m = nppcaUpdateM(model, expectations, varY, Y);
      model.mu = nppcaUpdateMu(model, expectations, varY, Y);
      
     case 'W'
      model.Cinv = nppcaUpdateCinv(model, expectations, varY, Y);
      model.W = nppcaUpdateW(model, expectations, varY, Y);
      
     case 'sigma'
      model.sigma = scg('nppcaSigmaObjective', model.sigma, options.optimiser, ...
                        'nppcaSigmaGradient', model, expectations, ...
                        varY, Y);
     
      
     
     case 'estep'
      model = nppcaRemoveRedundancy(model);
      expectations = nppcaEstep(model, expectations, varY, Y);  
    end
    if options.stepChecks
      [deltaL, oldL] = nppcaLikelihoodCheck(model, expectations, ...
                                            varY, Y, oldL, updateParam);
      maxDeltaL = max([maxDeltaL deltaL]);
    end
  end
  if ~options.stepChecks
    [maxDeltaL, oldL] = nppcaLikelihoodCheck(model, expectations, ...
                                             varY, Y, oldL, 'All params');
  end
  lastUpdate = order(end);
  order = randperm(length(params));
  while(lastUpdate == order(1))
    order = randperm(length(params));
  end
  fprintf('Iteration number: %d\n', counter);
 
  save tempResults model expectations counter
end

if counter >= options.maxIters
  fprintf('Warning maximum iterations exceeded.\n')
end

model = nppcaRemoveRedundancy(model);
expectations = nppcaEstep(model, expectations, varY, Y);  

save resultsCNS1 model expectations