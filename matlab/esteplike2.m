function f=esteplike2(mu,w,sigma,B,V,N,d,q,xmean,var,trvar)
%Likelihood Computes the likelihood of parameters for the E-step
%in the E-M algorythm for gMOSPPCAEM.
%w=k(1:d*q);
%sigma=k(q*d+1);
%mu=(k(q*d+2:end))';
x=reshape(w,d,q);
for i=1:N 

  z(i)=sum(log((sigma^2*ones(d,1)+B(:,i))));
  
  s(i)=((V(:,i)-mu)')*(((sigma^2*ones(d,1)+B(:,i)).^(-1)).*(V(:,i)-mu));
  s1(i)=(xmean(:,i)'*x'*(((sigma^2*ones(d,1)+B(:,i)).^(-1)).*((V(:,i)-mu))));
  for j=1:d
  s3(i,j)=trace(x(j,:)'*(((sigma^2+B(j,i))^(-1))*x(j,:))*var(:,:,i));          
  end
  s2=sum(s3,2)';
end
f=sum(z+trvar+s-2*s1+s2);
