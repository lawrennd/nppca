function f=mstep2(mu,w,sigma,B,V,N,d,q,xmean,var,trvar)
%MSTEP performs the M-step in the E_M algorythm of EMapprox;


%w=k(1:d*q);
%sigma=k(q*d+1);
%mu=(k(q*d+2:end))';
x=reshape(w,d,q);
%for i=1:N 
%  G(:,i)=((sigma^2*ones(d,1)+B(:,i)).^(-1)).*(V(:,i)-x*xmean(:,i));
  
%end
%Q=sum(((sigma^2*ones(d,N)+B).^(-1)),2);

%MU=(Q.\(sum(G,2)));
%MU=mu;

for i=1:N
 H(:,:,i)=(((sigma^2*ones(d,1)+B(:,i)).^(-1)).*(V(:,i)-mu))*xmean(:,i)';
  for j=1:d
   L(:,:,i,j)=((sigma^2+B(j,i)).^(-1))*var(:,:,i);
  end
end
h=sum(H,3);
l=sum(L,3);
for j=1:d
  NW(j,:)=h(j,:)/l(:,:,j);
end
f=NW(:)';
%for i=1:N
%  s(i)=((V(:,i)-mu')'*(((B(:,i)).^(-2)).*(V(:,i)-mu')));
%  s1(i)=(2*xmean(:,i)'*x'*(((B(:,i)).^(-2)).*(V(:,i)-mu')));
%  s2(i)=trace(diag((B(:,i)).^(-2))*(x*var(:,:,i)*x'));
%  z(i)=sum(B(:,i));
%end
%si=-(sum(s-s1+s2-z)/(N*d-2/sum(z)));
%if si>0
%SI=sqrt(si);   
%f=SI;
%else
%    fprintf('You are driving sigma to zero!\n');
%  f=sigma;
%end


