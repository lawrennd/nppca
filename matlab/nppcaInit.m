function model = nppcaInit(X, varX, latentDim);

% NPPCAINIT Initialise a Noisy probabilistic PCA model.

% NPPCA

G = cov(X);
[sigma2, V, s] = ppca(G,latentDim);
val = sqrt(sigma2)-mean(mean(varX));
if val < 0
  model.sigma = 1e-6;
else
  model.sigma = val;
end
model.W = V*diag(sqrt(s-model.sigma^2));  

model.mu = mean(X);
