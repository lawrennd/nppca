function W = nppcaUpdateW(model, expectations, B, V)

% NPPCAUPDATEW Update the W matrice for a noisy probabilistic PCA model.

% NPPCA

numData = size(V, 2);
d = size(V, 1);
q = size(expectations.x, 1);

for i=1:numData
 H(:,:,i)=(((model.sigma^2*ones(d,1)+B(:,i)).^(-1)).*(V(:,i)-model.mu))*expectations.x(:,i)';
  for j=1:d
   L(:,:,i,j)=((model.sigma^2+B(j,i)).^(-1))*expectations.xxT(:,:,i);
  end
end
h=sum(H,3);
l=sum(L,3);
for j=1:d
  W(j,:)=h(j,:)/l(:,:,j);
end


