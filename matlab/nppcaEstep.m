function expectations = nppcaEstep(model, expectations, B, X)

% NPPCAESTEP Estep for the noisy probabilistic PCA.

dataDim = size(X, 2);
latentDim = size(expectations.x, 1);
numData = size(X, 1);

for i=1:numData   %computes the inverses of the matrices M_n  
  Binv = diag((model.sigma^2*ones(dataDim,1) + B(i, :)').^(-1));
  Sigma_x = inv(model.W'*Binv*model.W + eye(latentDim));
  expectations.x(:,:,i) = Sigma_x*model.W'*Binv*(X(i, :) - model.mu)';
  expectations.xxT(:,:,i) = Sigma_x + expectations.x(:,:,i)*expectations.x(:,:,i)';
  expectations.xTx(i) = trace(expectations.xxT(:,:,i));
  expectations.logDetCov(i) = logdet(Sigma_x);
end 

