function m = nppcaUpdateM(model, expectations, B, X)

% NPPCAUPDATEM Update the latent mean for the noisy PPCA model.

% NPPCA

m = mean(expectations.x);

