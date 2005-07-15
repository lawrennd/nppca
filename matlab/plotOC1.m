%PLOTOC1 plots our OC1 data
% Set up the options.
options = nppcaOptions;
options.tol=1e-10;
options.maxIters=1e4;
options.display=2;
% Load first 20 points from OC1 data.
[probes, annotations, Y, varY] = nppcaLoadData('OC1');

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
  cModel = model.W*pdinv(model.Cinv)*model.W' + model.sigma* ...
           model.sigma*eye(dataDim);
  
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