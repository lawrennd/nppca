function mu = nppcaUpdateMu(model, expectations, B, X)

% NPPCAUPDATEMU Update the mean for the noisy PPCA model.

% NPPCA

numData = size(X, 1);
dataDim = size(X, 2);

for i=1:numData 
  G(:, i)=((model.sigma^2*ones(dataDim, 1)+B(i, :)').^(-1)).*(X(i, :)'-model.W*expectations.x(:, i));
end
Q=sum((model.sigma^2*ones(dataDim, numData)+B').^(-1), 2);

mu = Q.\sum(G, 2);
mu = mu';

