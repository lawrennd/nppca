function f = nppcaSigmaObjective(sigma, model, expectations, B, V);

% NPPCASIGMAOBJECTIVE Wrapper function for objective as a function of Sigma.

% NPPCA

model.sigma = sigma;
f = nppcaLikelihoodBound(model, expectations, B, V);