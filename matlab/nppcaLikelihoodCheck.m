function [deltaL, oldL] = nppcaLikelihoodCheck(model, expectations, varY, ...
                                       Y, oldL, param);

% NPPCALIKELIHOODCHECK Compute the difference in likelhoods.

L = nppcaLikelihoodBound(model, expectations, varY, Y);
deltaL = oldL - L;
oldL = L;
if deltaL < 0
  warning(['Likelihood drop of ' num2str(deltaL) ' after update of ' param '.']);
end
