% DEMNPPCA1 Simple demo of probabilistic PCA with noise.

% Fix a seed so that results are repeatable.
randn('seed', 2e5);
rand('seed', 2e5);

stepChecks = 0; % Whether to check likelihood after each update.
display = 1; % whether or not to display progress in a figure.
[probes, alpha, annotation, geneName] = gmosReadTxt(['../../gMOS/' ...
                    'data/signalOC1A.txt']);

alpha = alpha(1:100, :);
alphaLim = 1e-6;

alpha(find(alpha<1e-6)) = 1e-6;

Y = digamma(alpha);
varY = trigamma(alpha);

numData = size(Y, 1);
dataDim = size(Y, 2);
latentDim = 3;


% Number of EM iterations.
maxIters = 1000;
options = foptions; % optimisation options for Sigma.

% Initialise the model.
[model, expectations] = nppcaInit(Y, varY, latentDim);



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

% Compute the starting likelihood.
oldL = nppcaLikelihoodBound(model, expectations, varY, Y);
maxDeltaL = 1;
checkEvery = 10;
initIters = 10;
counter=0;

while  (maxDeltaL > 1e-5*checkEvery & counter < maxIters) ...
      | counter < initIters
  maxDeltaL = 0;
  counter=counter+1;
  if counter > initIters
    model.mu = nppcaUpdateMu(model, expectations, varY, Y);
    
    fprintf('Hmph ... mu is %2.4f %2.4f\n', model.mu(1), model.mu(2));
    if stepChecks
      [deltaL, oldL] = nppcaLikelihoodCheck(model, expectations, ...
                                            varY, Y, oldL, 'mu');
      maxDeltaL = max([maxDeltaL deltaL]);
    end
  end
  model.W = nppcaUpdateW(model, expectations, varY, Y);
  fprintf('Mmmmm ... W is %2.4f %2.4f\n', model.W(1), model.W(2));
  if stepChecks
    [deltaL, oldL] = nppcaLikelihoodCheck(model, expectations, ...
                                          varY, Y, oldL, 'W');
    maxDeltaL = max([maxDeltaL deltaL]);
  end
  %/~
  %  model.sigma = quasinew('nppcaSigmaObjective',model.sigma, options, ...
  %                         'nppcaSigmaGradient',model, expectations, varY, Y);
  %~/
  if counter < initIters | ~rem(counter, checkEvery)
    model.sigma = nppcaNewtonUpdateSigma(model,expectations,varY, Y);
    fprintf('Ahhh ... sigma is %2.4f\n', model.sigma);
    if stepChecks
      [deltaL, oldL] = nppcaLikelihoodCheck(model, expectations, ...
                                            varY, Y, oldL, 'sigma');
      maxDeltaL = max([maxDeltaL deltaL]);
    end
  end
  expectations = nppcaEstep(model, expectations, varY, Y);  
  if stepChecks
    [deltaL, oldL] = nppcaLikelihoodCheck(model, expectations, ...
                                  varY, Y, oldL, 'estep');
    maxDeltaL = max([maxDeltaL deltaL]);
  else
    if ~rem(counter, checkEvery)
      [maxDeltaL, oldL] = nppcaLikelihoodCheck(model, expectations, ...
                                               varY, Y, oldL, 'All params');
    else
      maxDeltaL = 1;
    end
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
  

