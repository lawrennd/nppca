function f=estepgradlike2(sigma, mu,w,B,V,N,d,q,xmean,var,trvar)
     %ESTEPGRADLIKE computes the gradient of the complete likelihood of estep;


x=reshape(w,d,q);
for i=1:N 
  F(:,:,i)=diag((sigma^2*ones(d,1)+B(:,i)).^(-1));
  Q(:,:,i)=(F(:,:,i).^2);
  s(i)=-2*((V(:,i)-mu)'*Q(:,:,i)*(V(:,i)-mu))*sigma;
  z(i)=2*trace((F(:,:,i)))*sigma;
  s2(i)=-2*trace(x'*Q(:,:,i)*x*var(:,:,i))*sigma;
  s1(i)=(4*xmean(:,i)'*x'*Q(:,:,i)*(V(:,i)-mu))*sigma;
  G(:,:,i)=-(F(:,:,i)*(2*(V(:,i)-mu)*xmean(:,i)'-2*x*var(:,:,i)));
  Sigma(i)=s(i)+z(i)+s1(i)+s2(i);
  Mu(:,i)=F(:,:,i)*(2*mu-2*V(:,i)+2*x*xmean(:,i));
end
GI=sum(G,3);
SI=sum(Sigma);
MU=sum(Mu,2);
f=SI;



