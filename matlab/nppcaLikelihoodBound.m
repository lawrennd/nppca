function f=nppcaLikelihoodBound(model, expectations, B, X)

% NPPCALIKELIHOODBOUND Gives the variational bound on the log likelihood for noisy PPCA.

% NPPCA

numData = size(X, 2);
dataDim = size(X, 1);

for i=1:numData   
  z(i)=sum(log((model.sigma^2*ones(dataDim,1)+B(:,i))));  
  s(i)=((X(:,i)-model.mu)')*(((model.sigma^2*ones(dataDim,1)+B(:,i)).^(-1)).*(X(:,i)-model.mu));
  s1(i)=(expectations.x(:,i)'*model.W'*(((model.sigma^2*ones(dataDim,1)+B(:,i)).^(-1)).*((X(:,i)-model.mu))));
  for j=1:dataDim
    s3(i,j)=trace(model.W(j,:)'*(((model.sigma^2+B(j,i))^(-1))*model.W(j,:))*expectations.xxT(:,:,i));          
  end
  s2=sum(s3,2)';
end
f=sum(z+expectations.xTx+s-2*s1+s2);
