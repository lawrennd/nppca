%EMapprox Given the output of a gMOS reading, 
%performs a PPCA-like maximum likelihood optimisation on a latent variable 
%model using an approximate E-M algorythm.It expects to read from a file 
%a matrix of readings, V, and a matrix of variances, B (readings along columns).

R=input('Which file should I read?','s');
load(R);
%Reads the file where the (column) vectors of means and variances are stored.
g=size(V);
N=g(2);  %N is the number of readings, usually ~10.
d=g(1);  %d is the dimensionality of the readings, usually ~10000
h=min(B(:));
B=B;
V=V;
model.mu=mean(V'); %initialises the mean vector to the empirical mean
             %Initialises the parameters of the LVM.
q=input('How many directions do you want to retain?\n');

tic;
  model.sigma=0.1;
  
  options = foptions;
  %options(9) = 1;
  options(1) = 1;
  options(2)=1e-4;
  options(3)=1e-4;
  options(14)=5000;
  %Initialises and starts the optimisation
  G=cov((V-model.mu'*ones(1,N))'./N);
  
  [s,z]=pca(G,q);
  %I use pca to initialise the factors;
  model.W=z(:)';  
  L=1;
  L1=0;   
  k=[model.W,model.sigma,model.mu]; %defines long vector with all parameters.
  t=0;
  x=reshape(k(1:d*q),d,q);
  model.sigma=k(d*q+1);
  model.mu=k(d*q+2:end);
  U=zeros(d,N);
  for i=1:N   %computes the inverses of the matrices M_n
    
    F(:,:,i)=diag((model.sigma^2*ones(d,1)+B(:,i)).^(-1));
    M(:,:,i)=inv(x'*F(:,:,i)*x+eye(q));
    expectations.x(:,:,i)=M(:,:,i)*x'*F(:,:,i)*(V(:,i)-model.mu');
    expectations.xxT(:,:,i)=M(:,:,i)+expectations.x(:,:,i)*expectations.x(:,:,i)';
    expectations.xTx(i)=trace(expectations.xxT(:,:,i));
  end 
  while abs((L1-L)/(abs(L)))>0.0001,
    t=t+1;
    fprintf('Step number');
    t
    
    %E-step
    L=feval('esteplike',k,B,V,N,d,q,expectations.x,expectations.xxT,expectations.xTx);
    %M-step
    k=feval('mstep',k,B,V,N,d,q,expectations.x,expectations.xxT,expectations.xTx);
    L1=feval('esteplike',k,B,V,N,d,q,expectations.x,expectations.xxT,expectations.xTx);
    if L1>L
      beep
      error('Che cazzo fai?\n');
         
    else
    end
    x=reshape(k(1:d*q),d,q);
    model.sigma=k(d*q+1);
    model.mu=k(d*q+2:end);
    for i=1:N %computes the inverses of the matrices M_n
      
      M(:,:,i)=inv(x'*F(:,:,i)*x+eye(q));
      expectations.x(:,i)=M(:,:,i)*x'*F(:,:,i)*(V(:,i)-model.mu');
      expectations.xxT(:,:,i)=M(:,:,i)+expectations.x(:,:,i)*expectations.x(:,:,i)';
      expectations.xTx(i)=trace(expectations.xxT(:,:,i));
      
      
    end 
   
  end
  model.W=k(1:d*q);
  model.sigma=k(d*q+1);
  model.mu=k(d*q+2:end);
  wML=reshape(model.W,d,q);
toc
beep
save results wML model.sigma model.mu
fprintf('Number of iterations used by EM algorythm\n')
t
fprintf('Output saved in results.mat\n')
