function g = nppcaSigmaCurvature(model, expectations, B, X)

% NPPCASIGMACURVATURE Wrapper function for curvature with respect to sigma.

numData = size(X, 1);
dataDim = size(X, 2);

for i = 1:numData 
  F = diag((model.sigma^2*ones(dataDim,1)+B(i,:)').^(-1));
  F2 = F.*F;
  F3 = F.*F2;
  s(i) = 8*((X(i, :) - model.mu)*F3*(X(i, :) - model.mu)')*model.sigma;
  z(i) = -4*trace(F2)*model.sigma;
  s2(i) = 8*trace(model.W'*F3*model.W*expectations.xxT(:,:,i))*model.sigma;
  s1(i) = (-16*expectations.x(i, :)*model.W'*F3*(X(i, :) - model.mu)')*model.sigma;
  Sigma(i) = s(i)+z(i)+s1(i)+s2(i);
end
g = sum(Sigma)*model.sigma + nppcaSigmaGradient(model, expectatiions, ...
                                               B, X)/model.sigma;