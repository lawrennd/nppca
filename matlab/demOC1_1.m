% DEMNPPCA1 Simple demo of probabilistic PCA with noise.

% Fix a seed so that results are repeatable.
randn('seed', 2e5);
rand('seed', 2e5);

stepChecks = 0; % Whether to check likelihood after each update.
display = 1; % whether or not to display progress in a figure.
[probes, alpha, annotation, geneName] = gmosReadTxt(['../../gMOS/' ...
                    'data/signalOC1A.txt']);

alpha = alpha(1:1000, :);
alphaLim = 1e-6;

alpha(find(alpha<1e-6)) = 1e-6;

Y = digamma(alpha(:, 1:2));
varY = trigamma(alpha(:, 1:2));

numData = size(Y, 1);
dataDim = size(Y, 2);
latentDim = 1;


% Number of EM iterations.
maxIters = 50;
options = foptions; % optimisation options for Sigma.

% Initialise the model.
model = nppcaInit(Y, varY, latentDim);

% Initialise the expectations
expectations.x = zeros(numData, latentDim);
expectations.xxT = zeros(latentDim, latentDim, numData);
expectations.xTx = zeros(numData, 1);


% Plot the data and ellipses indicating uncertainty.

if display
  plot(Y(:, 1), Y(:, 2), 'r.');
  hold on;
  theta = 0:0.02:2*pi;
  for i = 1:numData
    plot(Y(i, 1)+(sqrt(varY(i, 1)))*cos(theta),Y(i, 2)+(sqrt(varY(i, 2)))* ...
         sin(theta), 'r-');
  end
  
  % Plot the initialisation
  ppcaCentrePoint = plot(model.mu(1), model.mu(2), 'b+');
  cModel = model.W*model.W' + model.sigma*model.sigma*eye(dataDim);
  cModel = cModel(1:2, :);
  cModel = cModel(:, 1:2);
  [r, s] = eig(cModel);
  ellipse = r*[sqrt(s(1, 1))*cos(theta); sqrt(s(2, 2))*sin(theta)];
  ppcaCovHandle = line(model.mu(1)*ones(1,size(theta, 2))+ellipse(1, :), ...
                       model.mu(2)*ones(1,size(theta, 2))+ellipse(2, :));

end
% Do an E-step
expectations = nppcaEstep(model, expectations, varY, Y);

% Compute the starting likelihood.
oldL = nppcaLikelihoodBound(model, expectations, varY, Y);
maxDeltaL = 1;

counter=0;
while  maxDeltaL > 1e-5 & counter < maxIters
  maxDeltaL = 0;
  counter=counter+1;
    
  model.mu = nppcaUpdateMu(model, expectations, varY, Y);
  if stepChecks
    [deltaL, oldL] = nppcaLikelihoodCheck(model, expectations, ...
                                  varY, Y, oldL, 'mu');
    maxDeltaL = max([maxDeltaL deltaL]);
  end

  model.W = nppcaUpdateW(model, expectations, varY, Y);
  if stepChecks
    [deltaL, oldL] = nppcaLikelihoodCheck(model, expectations, ...
                                  varY, Y, oldL, 'W');
    maxDeltaL = max([maxDeltaL deltaL]);
  end
  
  model.sigma = nppcaNewtonUpdateSigma(model,expectations,varY, Y);
  if stepChecks
    [deltaL, oldL] = nppcaLikelihoodCheck(model, expectations, ...
                                  varY, Y, oldL, 'sigma');
    maxDeltaL = max([maxDeltaL deltaL]);
  end

  expectations = nppcaEstep(model, expectations, varY, Y);  
  if stepChecks
    [deltaL, oldL] = nppcaLikelihoodCheck(model, expectations, ...
                                  varY, Y, oldL, 'estep');
    maxDeltaL = max([maxDeltaL deltaL]);
  else
    [maxDeltaL, oldL] = nppcaLikelihoodCheck(model, expectations, ...
                                  varY, Y, oldL, 'estep');
  end
  
  
  % plots the data points with the respective (averaged) errors
  if display
    set(ppcaCentrePoint, 'Xdata', model.mu(1), 'Ydata', model.mu(2));
    cModel = model.W*model.W' + model.sigma*model.sigma*eye(dataDim);
    cModel = cModel(1:2, :);
    cModel = cModel(:, 1:2);
    [r, s] = eig(cModel);
    ellipse = r*[sqrt(s(1, 1))*cos(theta); sqrt(s(2, 2))*sin(theta)];
    set(ppcaCovHandle, 'Xdata', model.mu(1)*ones(1,size(theta, 2))+ellipse(1,:), 'Ydata', model.mu(2)*ones(1,size(theta, 2))+ellipse(2,:));
    drawnow
  end
  fprintf('Iteration number: %d\n', counter);
end
if counter >= maxIters
  fprintf('Warning maximum iterations exceeded.\n')
end
  

