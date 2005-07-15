function [model, expectations] = nppcaInit100(Y, varY, latentDim);

% NPPCAINIT Initialise a Noisy probabilistic PCA model.

% NPPCA
randn('seed', 2e5);
rand('seed', 2e5);
colordef white
% We request four latent dimensions, but model only uses 3.
latentDim = 7;

% Set up the options.
options = nppcaOptions;

Y(find(Y<1e-6)) = 1e-6;
Y=Y(1:100,:);
varY=varY(1:100,:);


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
%      model.sigma = nppcaNewtonUpdateLogSigma(model,expectations,varY, Y);
      model.sigma = scg('nppcaSigmaObjective', model.sigma, options.optimiser, ...
                        'nppcaSigmaGradient', model, expectations, varY, Y);
     
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
end

model = nppcaRemoveRedundancy(model);
expectations = nppcaEstep(model, expectations, varY, Y);  

save initialisarionCNS model expectations