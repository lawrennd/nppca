function [model, expectations] = nppcaInit(Y, varY, latentDim);

% NPPCAINIT Initialise a Noisy probabilistic PCA model.

% NPPCA

dataDim = size(Y, 2);
numData = size(Y, 1);
precYPlusSigma = 1./(varY); 
A = Y.*precYPlusSigma;
for i = 1:size(A, 1);
  A(i, :) = A(i, :)./sum(precYPlusSigma, 1);
end
model.mu = sum(A, 1);
for i = 1:size(Y, 1);
  V(i, :) = ((Y(i, :) - model.mu).^2.*precYPlusSigma(i, :))./ ...
            sum(precYPlusSigma,1);
end
model.sigma = mean(sum(V, 1)./sum(precYPlusSigma, 1));
model.W = randn(dataDim, latentDim)*0.1;
model.m = zeros(1, latentDim);
model.Cinv = eye(latentDim);
% Initialise the expectations
expectations.x = zeros(numData, latentDim);
expectations.xxT = zeros(latentDim, latentDim, numData);
expectations.xTx = zeros(numData, 1);
expectations = nppcaEstep(model, expectations, varY, Y);

