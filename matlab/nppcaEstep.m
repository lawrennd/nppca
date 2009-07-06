function expectations = nppcaEstep(model, expectations, varY, Y)

% NPPCAESTEP Estep for the noisy probabilistic PCA.

% NPPCA

dataDim = size(Y, 2);
latentDim = size(expectations.x, 2);
numData = size(Y, 1);

for i=1:numData   %computes the inverses of the matrices M_n  
  Binv = diag((model.sigma^2*ones(dataDim,1) + varY(i, :)').^(-1));
  % This inverse should really be done across a q dimensional matrix.
  SigmaInv = model.W'*Binv*model.W + model.Cinv;
  [Sigma_x, U] = pdinv(SigmaInv);
  expectations.logDetCov(i) = -logdet(SigmaInv, U);
  expectations.x(i, :) = (Sigma_x*(model.W'*Binv*(Y(i, :) - model.mu)'+model.Cinv*model.m'))';
  expectations.xxT(:,:,i) = Sigma_x ...
      + expectations.x(i, :)'*expectations.x(i, :);
%  expectations.xTx(i) = trace(expectations.xxT(:,:,i));
end 

