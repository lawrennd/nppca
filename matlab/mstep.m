function f=mstep(k,B,V,N,d,q,expectations.x,expectations.xxT,expectations.xTx)
%MSTEP performs the M-step in the E_M algorythm of EMapprox;


model.W=k(1:d*q);
model.sigma=k(q*d+1);
model.mu=(k(q*d+2:end))';
x=reshape(model.W,d,q);
for i=1:N 
  G(:,i)=((model.sigma^2*ones(d,1)+B(:,i)).^(-1)).*(V(:,i)-x*expectations.x(:,i));
  
end
Q=sum(((model.sigma^2*ones(d,N)+B).^(-1)),2);

MU=(Q.\(sum(G,2)));
%MU=model.mu;

for i=1:N
  H(:,:,i)=(((model.sigma^2*ones(d,1)+B(:,i)).^(-1)).*(V(:,i)-MU))*expectations.x(:,i)';
  for j=1:d
    L(:,:,i,j)=((model.sigma^2+B(j,i)).^(-1))*expectations.xxT(:,:,i);
  end
end
h=sum(H,3);
l=sum(L,3);
for j=1:d
  NW(j,:)=h(j,:)/l(:,:,j);
end

for i=1:N
  s(i)=((V(:,i)-MU)'*(V(:,i)-MU));
  s1(i)=(2*expectations.x(:,i)'*x'*(V(:,i)-MU));
  s2(i)=trace(x*expectations.xxT(:,:,i)*x');
  z(i)=sum(B(:,i));
end
si=(sum(s-s1+s2+z)/(N*d));
if si>0
SI=model.sigma;  
  f=[NW(:)',SI,MU'];
else
  fprintf('You are driving model.sigma to zero!\n');
  f=[NW(:)',model.sigma,MU'];
end


