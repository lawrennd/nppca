function f = newtonUpdate(model, expectations, B, X)
%NEWTONUPDATE performs a Newton update for the sigma dependence of
%the likelihood.
f = model.sigma - newtonStep(model, expectations, B, X);