function g = nppcaSigmaGradient(sigma, model, expectations, B, X)

% NPPCASIGMAGRADIENT Wrapper function for gradient with respect to sigma.

numData = size(X, 1);
dataDim = size(X, 2);

for i = 1:numData 
  F = diag((model.sigma^2*ones(dataDim,1)+B(i,:)').^(-1));
  F2 = F.*F;
  s(i) = -2*((X(i, :) - model.mu)*F2*(X(i, :) - model.mu)')*model.sigma;
  z(i) = 2*trace(F)*model.sigma;
  s2(i) = -2*trace(model.W'*F2*model.W*expectations.xxT(:,:,i))*model.sigma;
  s1(i) = (4*expectations.x(i, :)*model.W'*F2*(X(i, :) - model.mu)')*model.sigma;
  Sigma(i) = s(i)+z(i)+s1(i)+s2(i);
end
g = 0.5*sum(Sigma);



