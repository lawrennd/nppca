%PLOTOC1MANC plots the manky dataset
% Load points from OC1 data.

Y=load('../../gMOS/data/MU_setB_e_mmgmos.txt');
varY=load('../../gMOS/data/MU_setB_large_var_mmgmos.txt');

% Initialise the model --- reset to PCA.
[model, expectations] = nppcaInit(Y, varY, latentDim);

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

 
  fprintf('Press any key to continue\n');
  pause