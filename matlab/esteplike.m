function f=esteplike(k,B,V,N,d,q,expectations.x,expectations.xxT,expectations.xTx)
%Likelihood Computes the likelihood of parameters for the E-step
%in the E-M algorythm for gMOSPPCAEM.
model.W=k(1:d*q);
model.sigma=k(q*d+1);
model.mu=(k(q*d+2:end))';
x=reshape(model.W,d,q);
for i=1:N 

  z(i)=sum(log((model.sigma^2*ones(d,1)+B(:,i))));
  
  s(i)=((V(:,i)-model.mu)')*(((model.sigma^2*ones(d,1)+B(:,i)).^(-1)).*(V(:,i)-model.mu));
  s1(i)=trace(expectations.x(:,i)'*x'*(((model.sigma^2*ones(d,1)+B(:,i)).^(-1)).*((V(:,i)-model.mu))));
  for j=1:d
  s3(i,j)=trace(x(j,:)'*(((model.sigma^2+B(j,i)).^(-1))*x(j,:))*expectations.xxT(:,:,i));          
  end
  s2=sum(s3,2)';
end
f=sum(z+expectations.xTx+s-2*s1+s2)+N*d*log(2*pi)/2;
