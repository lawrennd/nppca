% DEMNPPCA1 Simple demo of probabilistic PCA with noise.

% Fix a seed so that results are repeatable.
randn('seed', 2e5);
rand('seed', 2e5);

display = 0;
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
maxIters = 9;
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
while  maxDeltaL > 1e-3 & counter < maxIters
  maxDeltaL = 0;
  counter=counter+1;
    
  model.mu = nppcaUpdateMu(model, expectations, varY, Y);
  L = nppcaLikelihoodBound(model, expectations, varY, Y);
  deltaL = oldL - L;
  oldL = L;
  if deltaL < 0
    warning(['Likelihood drop of ' num2str(deltaL) ' after update of mu.']);
  end
  maxDeltaL = max([maxDeltaL deltaL]);

  model.W = nppcaUpdateW(model, expectations, varY, Y);
  L = nppcaLikelihoodBound(model, expectations, varY, Y);
  deltaL = oldL - L;
  oldL = L;
  if deltaL < 0
    warning(['Likelihood drop of ' num2str(deltaL) ' after update of w.']);
  end
  maxDeltaL = max([maxDeltaL deltaL]);

  model.sigma = nppcaNewtonUpdateSigma(model,expectations,varY, Y);
  L = nppcaLikelihoodBound(model, expectations, varY, Y);
  deltaL = oldL - L;
  oldL = L;
  if deltaL < 0
    warning(['Likelihood drop of ' num2str(deltaL) ' after update of sigma.']);
  end
  maxDeltaL = max([maxDeltaL deltaL]);

  expectations = nppcaEstep(model, expectations, varY, Y);  
  L = nppcaLikelihoodBound(model, expectations, varY, Y);
  deltaL = oldL - L;
  oldL = L;
  if deltaL < 0
    warning(['Likelihood drop of ' num2str(deltaL) ' after E step.']);
  end
  maxDeltaL = max([maxDeltaL deltaL]);
  
  
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
  

