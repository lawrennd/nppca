function [model, expectations] = nppcaPpcaInit(Y, varY, latentDim, numData)

G = cov(Y);
[sigma2, V, s] = ppca(G, latentDim);
model.sigma=sqrt(sigma2);
model.m = zeros(1, latentDim);
model.Cinv = eye(latentDim);
model.W = V*diag(sqrt(s-model.sigma^2));  
model.mu = mean(Y);
% Initialise the expectations
expectations.x = zeros(numData, latentDim);
expectations.xxT = zeros(latentDim, latentDim, numData);
expectations.xTx = zeros(numData, 1);
expectations = nppcaEstep(model, expectations, varY, Y);

