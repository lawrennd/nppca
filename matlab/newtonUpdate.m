function f = newtonUpdate(model, expectations, varY, Y)
%NEWTONUPDATE performs a Newton update for the sigma dependence of
%the likelihood.
f = model.sigma - nppcaSigmaGradient(model, expectations, varY, Y)/nppcaSigmaCurvature(model, expectations, varY, Y);