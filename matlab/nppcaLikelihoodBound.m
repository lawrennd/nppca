function f=nppcaLikelihoodBound(model, expectations, B, X)

% NPPCALIKELIHOODBOUND Gives the variational bound on the log likelihood for noisy PPCA.

% NPPCA

numData = size(X, 1);
dataDim = size(X, 2);

z = zeros(numData, 1);
s = zeros(numData, 1);
s1 = zeros(numData, 1);
s3 = zeros(numData, dataDim);

sigma2 = model.sigma*model.sigma;
for i = 1:numData   
  xHat = X(i, :) - model.mu;
  Binv = 1./(sigma2*ones(dataDim,1) + B(i, :)');
  z(i) = -sum(log(Binv));  
  s(i) = xHat*(Binv.*xHat');
  s1(i) = expectations.x(i, :)*model.W'*(Binv.*xHat');
  for j = 1:dataDim
    s3(i, j) = Binv(j)*model.W(j,:)*expectations.xxT(:,:,i)*model.W(j,:)';     
  end
  s2 = sum(s3,2);
end
f = sum(z + expectations.xTx + s - 2*s1 + s2);
