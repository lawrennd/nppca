function f=mstep(k,B,V,N,d,q,xmean,var,trvar)
%MSTEP performs the M-step in the E_M algorythm of EMapprox;


w=k(1:d*q);
sigma=k(q*d+1);
mu=(k(q*d+2:end))';
x=reshape(w,d,q);
for i=1:N 
  G(:,i)=((sigma^2*ones(d,1)+B(:,i)).^(-1)).*(V(:,i)-x*xmean(:,i)');
  
end
Q=sum(((sigma^2*ones(d,N)+B).^(-1)),2);

MU=(Q.\(sum(G,2)));
for i=1:N
  H(:,:,i)=(((sigma^2*ones(d,1)+B(:,i)).^(-1)).*(V(:,i)-MU))*xmean(:,i)';
  for j=1:d
    L(:,:,i,j)=((sigma^2+B(j,i)).^(-1))*var(:,:,i);
  end
end
h=sum(H,3);
l=sum(L,3);
for j=1:d
  NW(j,:)=h(j,:)/l(:,:,j);
end

for i=1:N
  s(i)=((V(:,i)-MU)'*(V(:,i)-MU));
  s1(i)=(2*xmean(:,i)'*NW'*(V(:,i)-MU));
  s2(i)=trace(NW*var(:,:,i)*NW');
  z(i)=sum(B(:,i));
end
si=(sum(s-s1+s2-z)/(N*d));
if si>0
SI=sqrt(si);  
  f=[NW(:)',SI,MU'];
else
  fprintf('You are driving sigma to zero!\n');
  f=[NW(:)',sigma,MU'];
end


