function f=esteplike2(model, expectations, B, V)

% ESTEPLIKE2 Computes negative log expectation of the joint distribution.

N = size(V, 2);
d = size(V, 1);
q = size(expectations.x, 1);
x=reshape(model.W,d,q);
for i=1:N 
  
  z(i)=sum(log((model.sigma^2*ones(d,1)+B(:,i))));  
  s(i)=((V(:,i)-model.mu)')*(((model.sigma^2*ones(d,1)+B(:,i)).^(-1)).*(V(:,i)-model.mu));
  s1(i)=(expectations.x(:,i)'*model.W'*(((model.sigma^2*ones(d,1)+B(:,i)).^(-1)).*((V(:,i)-model.mu))));
  for j=1:d
    s3(i,j)=trace(model.W(j,:)'*(((model.sigma^2+B(j,i))^(-1))*model.W(j,:))*expectations.xxT(:,:,i));          
  end
  s2=sum(s3,2)';
end
f=sum(z+expectations.xTx+s-2*s1+s2);
