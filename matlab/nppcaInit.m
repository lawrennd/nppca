function model = nppcaInit(X, latentDim);

% NPPCAINIT Initialise a Noisy probabilistic PCA model.

% NPPCA

G = cov(X);
[sigma2, z, s] = ppca(G,latentDim);
model.sigma = sqrt(sigma2);
model.W = sqrt(s-model.sigma^2)*z(:)';  
model.W = model.W';
model.mu = mean(X)';
