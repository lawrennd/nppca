function W = nppcaUpdateW(model, expectations, B, X)

% NPPCAUPDATEW Update the W matrice for a noisy probabilistic PCA model.

% NPPCA

numData = size(X, 1);
dataDim = size(X, 2);

for i = 1:numData
  H(:,:,i) = (((model.sigma^2*ones(dataDim, 1)+B(i, :)').^(-1))...
            .*(X(i, :) - model.mu)')*expectations.x(i, :);
  for j = 1:dataDim
    L(:,:,i,j) = ((model.sigma^2 + B(i, j)).^(-1))*expectations.xxT(:,:,i);
  end
end
h = sum(H,3);
l = sum(L,3);
for j = 1:dataDim
  W(j,:) = h(j,:)/l(:,:,j);
end


