function sigma = nppcaNewtonUpdateLogSigma(model, expectations, varY, Y)

% NPPCANEWTONUPDATELOGSIGMA performs a Newton update for the logsigma dependence of the likelihood.

% NPPCA
sigma = model.sigma;
deltaL1 = 1;
maxIters = 10;
counter2 = 0;
while deltaL1 > 1e-3 & counter2 <= 10
  deltaL = -1;
  counter2 = counter2 + 1;
  oldL = nppcaSigmaObjective(sigma, ...
                             model, expectations, varY, Y);
  step = newtonStep(sigma, model, expectations, varY, Y);
  counter = 0;
  while deltaL < 0
    counter = counter + 1;
    L = nppcaSigmaObjective(sigma-step, ...
                            model, expectations, varY, Y);
    deltaL = oldL - L;
    step = step/2;
  end
  deltaL1 = deltaL;
  sigma = sigma - 2*step;
  %/~
  %fprintf('Now Sigma is %2.4f\n', sigma);
  %~/
end
function step = newtonStep(sigma, model, expectations, varY, Y)

% NEWTONSTEP Wrapper function for gradient with respect to sigma.

if sigma < eps
  sigma = eps;
end
eta = 0.1;
numData = size(Y, 1);
dataDim = size(Y, 2);


for i = 1:numData 
  F = diag((sigma^2*ones(dataDim,1)+varY(i,:)').^(-1));
  F2 = F.*F;
  F3 = F2.*F;
  %computes the gradient with respect to sigma
  s(i) = -2*((Y(i, :) - model.mu)*F2*(Y(i, :) - model.mu)')*sigma;
  z(i) = 2*trace(F)*sigma;
  s2(i) = -2*trace(model.W'*F2*model.W*expectations.xxT(:,:,i))*sigma;
  s1(i) = (4*expectations.x(i, :)*model.W'*F2*(Y(i, :) - model.mu)')*sigma;
  Sigma(i) = s(i)+z(i)+s1(i)+s2(i);

  %computes the curvature with respect to sigma
  S(i) = 8*((Y(i, :) - model.mu)*F3*(Y(i, :) - model.mu)')*sigma;
  Z(i) = -4*trace(F2)*sigma;
  S2(i) = 8*trace(model.W'*F3*model.W*expectations.xxT(:,:,i))*sigma;
  S1(i) = (-16*expectations.x(i, :)*model.W'*F3*(Y(i, :) - model.mu)')*sigma;
  curvature(i) = S(i)+Z(i)+S1(i)+S2(i);
end
l=log(sigma);
Grad = 0.5*sum(Sigma);
%Grad = nppcaSigmaGradient(sigma, model, expectations, varY, Y)*sigma/2;
Curv = (0.5*sum(curvature)*sigma + Grad/sigma)*sigma*sigma/4;
%if Curv < 0
%  step = sign(grad)*model.sigma/2;
%  fprintf('Curvature is %2.4f. Gradient %f, setting step to %f\n', ...
%          Curv, Grad, step);
%end
%/~

%~/
step.l=Grad*exp(l)/(Curv*exp(2*l)+Grad*exp(l));

step=exp(l-step.l);