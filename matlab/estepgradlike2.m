function f=estepgradlike2(sigma, model, expectations, B, V)

% ESTEPGRADLIKE computes the gradient of the complete likelihood of estep;

N = size(V, 2);
d = size(V, 1);
q = size(expectations.x, 1);
x=reshape(model.W,d,q);
for i=1:N 
  F(:,:,i)=diag((model.sigma^2*ones(d,1)+B(:,i)).^(-1));
  Q(:,:,i)=(F(:,:,i).^2);
  s(i)=-2*((V(:,i)-model.mu)'*Q(:,:,i)*(V(:,i)-model.mu))*model.sigma;
  z(i)=2*trace((F(:,:,i)))*model.sigma;
  s2(i)=-2*trace(x'*Q(:,:,i)*x*expectations.xxT(:,:,i))*model.sigma;
  s1(i)=(4*expectations.x(:,i)'*x'*Q(:,:,i)*(V(:,i)-model.mu))*model.sigma;
  G(:,:,i)=-(F(:,:,i)*(2*(V(:,i)-model.mu)*expectations.x(:,i)'-2*x*expectations.xxT(:,:,i)));
  Sigma(i)=s(i)+z(i)+s1(i)+s2(i);
  Mu(:,i)=F(:,:,i)*(2*model.mu-2*V(:,i)+2*x*expectations.x(:,i));
end
GI=sum(G,3);
SI=sum(Sigma);
MU=sum(Mu,2);
f=SI;



