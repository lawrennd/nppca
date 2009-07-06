function [model, expectations] = nppcaFAInit(Y, varY, latentDim)

% NPPCAFAINIT Initialises NPPCA with factor analysis with constant errors

% NPPCA

numData=size(Y,1);
avgVarY=mean(varY,1);
scaledY = Y/diag(sqrt(avgVarY));%scale the rows to transform FA to PCA
G = cov(scaledY);
[sigma2, V, s] = ppca(G, latentDim);
% Initialise the model parameters
model.sigma=sqrt(sigma2);
model.m = zeros(1, latentDim);
model.Cinv = eye(latentDim);
model.W = diag(sqrt(avgVarY))*(V*diag(sqrt(s-model.sigma^2)));  
%initialise mean to weighted mean
weightedY=zeros(size(Y,1),size(Y,2));
for i=1:numData
    weightedY(i,:)=diag(varY(i,:))*Y(i,:)';
end
model.mu = sum(weightedY,1)/diag(sum(varY,1));
% Initialise the expectations
expectations.x = zeros(numData, latentDim);
expectations.xxT = zeros(latentDim, latentDim, numData);
expectations.xTx = zeros(numData, 1);
expectations = nppcaEstep(model, expectations, varY, Y);

