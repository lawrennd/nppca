function model = nppcaInit(X, latentDim);

% NPPCAINIT Initialise a Noisy probabilistic PCA model.

% NPPCA

G = cov(X);
[sigma2, V, s] = ppca(G,latentDim);
model.sigma = sqrt(sigma2);
model.W = V*diag(sqrt(s-model.sigma^2));  
model.mu = mean(X);
