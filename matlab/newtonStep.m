function g = newtonStep(model, expectations, B, X)

% NPPCASIGMAGRADIENT Wrapper function for gradient with respect to sigma.

numData = size(X, 1);
dataDim = size(X, 2);

for i = 1:numData 
  F = diag((model.sigma^2*ones(dataDim,1)+B(i,:)').^(-1));
  F2 = F.*F;
  F3 = F2.*F;
  %computes the gradient with respect to sigma
  s(i) = -2*((X(i, :) - model.mu)*F2*(X(i, :) - model.mu)')*model.sigma;
  z(i) = 2*trace(F)*model.sigma;
  s2(i) = -2*trace(model.W'*F2*model.W*expectations.xxT(:,:,i))*model.sigma;
  s1(i) = (4*expectations.x(i, :)*model.W'*F2*(X(i, :) - model.mu)')*model.sigma;
  Sigma(i) = s(i)+z(i)+s1(i)+s2(i);
  %computes the curvature with respect to sigma
  S(i) = 8*((X(i, :) - model.mu)*F3*(X(i, :) - model.mu)')*model.sigma;
  Z(i) = -4*trace(F2)*model.sigma;
  S2(i) = 8*trace(model.W'*F3*model.W*expectations.xxT(:,:,i))*model.sigma;
  S1(i) = (-16*expectations.x(i, :)*model.W'*F3*(X(i, :) - model.mu)')*model.sigma;
  SIgma(i) = S(i)+Z(i)+S1(i)+S2(i);
end
Grad = sum(Sigma);
Curv = sum(SIgma)*model.sigma + sum(Sigma)/model.sigma;

step=Grad/Curv;
if abs(step) > 0.1*model.sigma
  g = 0.2*step;
else
  g = Grad/Curv;
end
