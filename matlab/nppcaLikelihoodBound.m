function f=nppcaLikelihoodBound(model, expectations, varY, Y)

% NPPCALIKELIHOODBOUND Gives the variational bound on the log likelihood for noisy PPCA.

% NPPCA

numData = size(Y, 1);
dataDim = size(Y, 2);

z = zeros(numData, 1);
s = zeros(numData, 1);
s1 = zeros(numData, 1);
s3 = zeros(numData, dataDim);
trxHatCxHatT = zeros(numData, 1);
sigma2 = model.sigma*model.sigma;
for i = 1:numData   
  yHat = Y(i, :) - model.mu;
  Binv = 1./(sigma2 + varY(i, :)');
  z(i) = -sum(log(Binv));  
  s(i) = yHat*(Binv.*yHat');
  s1(i) = expectations.x(i, :)*model.W'*(Binv.*yHat');
  for j = 1:dataDim
    s3(i, j) = Binv(j)*model.W(j,:)*expectations.xxT(:,:,i)*model.W(j,:)';     
  end
  s2 = sum(s3,2);
  trxHatCxHatT(i) = traceProduct(expectations.xxT(:, :, i), model.Cinv) ...
      - 2*expectations.x(i, :)*model.Cinv*model.m' + model.m*model.Cinv*model.m';
end
f = 0.5*sum(z + trxHatCxHatT + s - 2*s1 + s2) + ...
    - 0.5*numData*logdet(model.Cinv);
f = f - 0.5*sum(expectations.logDetCov);

