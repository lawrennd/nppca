% DEMOC1Igata3manc Simple demo of probabilistic PCA with noise on
% reduced gata3 dataset with variance reduced by 4.

% Fix a seed so that results are repeatable.
randn('seed', 2e5);
rand('seed', 2e5);
colordef white
% We request four latent dimensions, but model only uses 3.
latentDim = 7;

% Set up the options.
options = nppcaOptions;

% Load first 20 points from OC1 data.
Y=load('../../gMOS/data/MU_setB_e_mmgmos.txt');
varY=load('../../gMOS/data/MU_setB_large_var_mmgmos.txt');

 YLim = 1e-6;
   
   Y(find(Y<1e-6)) = 1e-6;
Y=Y(5001:5500,:);
varY=varY(5001:5500,:);

% Initialise the model --- reset to PCA.
[model, expectations] = nppcaInit(Y, varY, latentDim);

if options.display > 1
  % Plot the data and ellipses indicating uncertainty.
  % We use only first two data dimensions.
  numData = size(Y, 1);
  dataDim = size(Y, 2);
  set(gca, 'fontsize', 20)
  set(gca, 'fontname', 'helvetica')
  
  % Plot points
  plot(Y(:, 1), Y(:, 2), 'r.');
  xlabel('y_1')
  ylabel('y_2')
  hold on;

  % Plot elipses
  theta = 0:0.02:2*pi;
  for i = 1:numData
    plot(Y(i, 1)+(sqrt(varY(i, 1)))*cos(theta),Y(i, 2)+(sqrt(varY(i, 2)))* ...
         sin(theta), 'r-');
  end
  origXlim = get(gca, 'xlim');
  origYlim = get(gca, 'ylim');
  line(origXlim, [0 0], 'color', [0 0 0])
  line([0 0], origYlim, 'color', [0 0 0])
  % Plot the initialisation
  mu = model.mu + model.m*model.W';
  ppcaCentrePoint = plot(mu(1), mu(2), 'b+');
  set(ppcaCentrePoint, 'lineWidth', 2, 'markersize', 10);
  cModel = model.W*pdinv(model.Cinv)*model.W' + model.sigma*model.sigma*eye(dataDim);
  cModel = cModel(1:2, :);
  cModel = cModel(:, 1:2);
  [r, s] = eig(cModel);
  
  ellipse = r*[sqrt(s(1, 1))*cos(theta); sqrt(s(2, 2))*sin(theta)];
  ppcaCovHandle = line(mu(1)*ones(1,size(theta, 2))+ellipse(1, :), ...
                       mu(2)*ones(1,size(theta, 2))+ellipse(2, :), ...
                       'linewidth', 2, 'color', [0 0 1]);
  print -depsc figure1
  fprintf('Press any key to continue\n');
  pause
end

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
      model.sigma = nppcaNewtonUpdateLogSigma(model,expectations,varY, Y);
%      model.sigma = scg('nppcaSigmaObjective', model.sigma, options.optimiser, ...
%                        'nppcaSigmaGradient', model, expectations, varY, Y);
     
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
  
  % plots the data points with the respective (averaged) errors
  if options.display > 1
    mu = model.m*model.W' + model.mu;
    set(ppcaCentrePoint, 'Xdata', mu(1), 'Ydata', mu(2));
    cModel = model.W*pdinv(model.Cinv)*model.W' + model.sigma*model.sigma*eye(dataDim);
    cModel = cModel(1:2, :);
    cModel = cModel(:, 1:2);
    [r, s] = eig(cModel);
    ellipse = r*[sqrt(s(1, 1))*cos(theta); sqrt(s(2, 2))*sin(theta)];
    set(ppcaCovHandle, 'Xdata', mu(1)*ones(1,size(theta, 2))+ellipse(1,:), 'Ydata', mu(2)*ones(1,size(theta, 2))+ellipse(2,:));
    drawnow
   
  end
  fprintf('Iteration number: %d\n', counter);
end


model = nppcaRemoveRedundancy(model);
expectations = nppcaEstep(model, expectations, varY, Y);  

save resultsGata3 model expectations Y varY
figure, nppcaProfilePlotter(model,expectations,Y, varY,484,'gata3')