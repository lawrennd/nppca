function f = nppcaSigmaObjective(sigma, model, expectations, B, V);

% NPPCASIGMAOBJECTIVE Wrapper function for bjective as a function of Sigma.

% NPPCA

model.sigma = sigma;
f = nppcaLikelihoodBound(model, expectations, B, V);