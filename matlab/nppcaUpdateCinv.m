function Cinv = nppcaUpdateCinv(model, expectations, B, X)

% NPPCAUPDATECINV Update the latent precision for the noisy PPCA model.

% NPPCA

C = mean(expectations.xxT, 3) ...
    - 2*mean(expectations.x)'*model.m ...
    + model.m'*model.m;

Cinv = pdinv(C);