function model = nppcaRemoveRedundancy(model)

% NPPCAREMOVEREDUNDANCY Remove the redundancy in m and Cinv.

% NPPCA

latentDim = size(model.W, 2);
model.mu = model.mu + model.m*model.W';
C = model.W*pdinv(model.Cinv)*model.W';
C = C/2 +C'/2;
[U, V] = eig(C);
V = diag(V);
[V, index] = sort(V);
index = index(end:-1:1);
V = V(end:-1:1);
U = U(:, index);

model.W = U(:, 1:latentDim)*diag(sqrt(V(1:latentDim)));
model.Cinv = eye(size(model.Cinv));
model.m = zeros(size(model.m));
      