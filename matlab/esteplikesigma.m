function f=esteplikesigma(sigma, model, expectations, B, V);

model.sigma = sigma;
f = nppcaLikelihoodBound(model, expectations, B, V);