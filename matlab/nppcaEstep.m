function expectations = nppcaEstep(model, expectations, B, logX)

% NPPCAESTEP Estep for the noisy probabilistic PCA.

dataDim = size(logX, 1);
latentDim = size(expectations.x, 1);
numData = size(logX, 2);
x=reshape(model.W,dataDim,latentDim);    
for i=1:numData   %computes the inverses of the matrices M_n  
  Binv = diag((model.sigma^2*ones(dataDim,1) + B(:,i)).^(-1));
  Sigma_x = inv(x'*Binv*x+eye(latentDim));
  expectations.x(:,:,i) = Sigma_x*x'*Binv*(logX(:,i)-model.mu);
  expectations.xxT(:,:,i) = Sigma_x + expectations.x(:,:,i)*expectations.x(:,:,i)';
  expectations.xTx(i) = trace(expectations.xxT(:,:,i));
  expectations.logDetCov(i) = logdet(Sigma_x);
end 

