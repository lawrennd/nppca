% DEMNPPCA1 Simple demo of probabilistic PCA with noise.

% Fix a seed so that results are repeatable.
randn('seed', 2e5);
rand('seed', 2e5);

% How many samples to take
numData = 10;
dataDim = 2; % data dimension.
latentDim = 1;

% Number of EM iterations.
maxIters = 10;
options = foptions; % optimisation options for Sigma.
options(9) = 1;
options(1) = 1;

logVals = zeros(numData, dataDim);
vals = gammarnd(3, 0.3, numData, 1);
logVals(:, 1) = digamma(vals);
logVals(:, 2) = logVals(:, 1)*2;
varVals = trigamma([vals vals.*vals]);
logVals = logVals + randn(size(logVals)).*sqrt(varVals);
X = logVals;
B = varVals;
B = B/10;

% Initialise the model.
model = nppcaInit(X, latentDim);

% Initialise the expectations
expectations.x = zeros(latentDim, 1, numData);
expectations.xxT = zeros(latentDim, latentDim, numData);
expectations.xTx = zeros(1, numData);


% Plot the data and ellispses indicating uncertainty.
plot(X(:, 1), X(:, 2), 'r.');
hold on;
theta = 0:0.02:2*pi;
for i = 1:numData
  plot(X(i, 1)+(sqrt(B(i, 1)))*cos(theta),X(i, 2)+(sqrt(B(i, 2)))* ...
       sin(theta), 'r-');
end

% Plot the initialisation
ppcaCentrePoint = plot(model.mu(1), model.mu(2), 'b+');
z = model.W/norm(model.W);
s = eig(model.W*model.W');
r = [z [-z(2); z(1)]];
ellipse = r*[sqrt(norm(model.W))*cos(theta); model.sigma*sin(theta)];
ppcaCovHandle = line(model.mu(1)*ones(1,size(theta, 2))+ellipse(1, :),model.mu(2)*ones(1,size(theta, 2))+ellipse(2,:));

% Do an E-step
expectations = nppcaEstep(model, expectations, B, X);


% Compute the starting likelihood.
oldL = nppcaLikelihoodBound(model, expectations, B, X);
deltaL = 1;

counter=0;
while  deltaL > 1e-5 & counter < maxIters
  counter=counter+1;
  
  model.mu = nppcaUpdateMu(model, expectations, B, X);
  L = nppcaLikelihoodBound(model, expectations, B, X);
  deltaL = oldL - L;
  oldL = L;
  if deltaL < 0
    warning(['Likelihood drop of ' num2str(deltaL) ' after update of mu.']);
  end

  model.W = nppcaUpdateW(model, expectations, B, X);
  L = nppcaLikelihoodBound(model, expectations, B, X);
  deltaL = oldL - L;
  oldL = L;
  if deltaL < 0
    warning(['Likelihood drop of ' num2str(deltaL) ' after update of w.']);
  end
    
  model.sigma = scg('nppcaSigmaObjective',model.sigma, options, ...
                  'nppcaSigmaGradient',model, expectations, B, X);
  L = nppcaLikelihoodBound(model, expectations, B, X);
  deltaL = oldL - L;
  oldL = L;
  if deltaL < 0
    warning(['Likelihood drop of ' num2str(deltaL) ' after update of sigma.']);
  end
  
  expectations = nppcaEstep(model, expectations, B, X);  
  L = nppcaLikelihoodBound(model, expectations, B, X);
  deltaL = 1;
  oldL = L;
  
  
  % plots the data points with the respective (averaged) errors
  z = model.W/norm(model.W);
  s = eig(model.W'*model.W);
  
  set(ppcaCentrePoint, 'Xdata', model.mu(1), 'Ydata', model.mu(2));
  r = [z, [-z(2); z(1)]];
  ellipse = r*[sqrt(s)*cos(theta); model.sigma*sin(theta)];
  set(ppcaCovHandle, 'Xdata', model.mu(1)*ones(1,size(theta, 2))+ellipse(1,:), 'Ydata', model.mu(2)*ones(1,size(theta, 2))+ellipse(2,:));
  fprintf('Iteration number: %d\n', counter);
  drawnow
end
if counter >= maxIters
  fprintf('Warning maximum iterations exceeded.')
end
save results model
  

