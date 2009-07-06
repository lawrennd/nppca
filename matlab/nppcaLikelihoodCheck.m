function [deltaL, L] = nppcaLikelihoodCheck(model, expectations, varY, ...
                                       Y, oldL, param);

% NPPCALIKELIHOODCHECK Compute the difference in likelhoods.

% NPPCA

L = nppcaLikelihoodBound(model, expectations, varY, Y);
deltaL = oldL - L;
oldL = L;
fprintf('Likelihood change of with update of %s: %2.4f\n', param, deltaL)
if deltaL < 0
  warning(['Likelihood drop of ' num2str(deltaL) ' after update of ' param '.']);
end
