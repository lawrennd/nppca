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
mu=mean(V'); %initialises the mean vector to the empirical mean
             %Initialises the parameters of the LVM.
q=input('How many directions do you want to retain?\n');

tic;
  sigma=0.1;
  
  options = foptions;
  %options(9) = 1;
  options(1) = 1;
  options(2)=1e-4;
  options(3)=1e-4;
  options(14)=5000;
  %Initialises and starts the optimisation
  G=cov((V-mu'*ones(1,N))'./N);
  
  [s,z]=pca(G,q);
  %I use pca to initialise the factors;
  w=z(:)';  
  L=1;
  L1=0;   
  k=[w,sigma,mu]; %defines long vector with all parameters.
  t=0;
  x=reshape(k(1:d*q),d,q);
  sigma=k(d*q+1);
  mu=k(d*q+2:end);
  U=zeros(d,N);
  for i=1:N   %computes the inverses of the matrices M_n
    
    F(:,:,i)=diag((sigma^2*ones(d,1)+B(:,i)).^(-1));
    M(:,:,i)=inv(x'*F(:,:,i)*x+eye(q));
    xmean(:,:,i)=M(:,:,i)*x'*F(:,:,i)*(V(:,i)-mu');
    var(:,:,i)=M(:,:,i)+xmean(:,:,i)*xmean(:,:,i)';
    trvar(i)=trace(var(:,:,i));
  end 
  while abs(L1-L)>0.0001*sqrt(d),
    t=t+1;
    fprintf('Step number');
    t
    
    %E-step
    L=feval('esteplike',k,B,V,N,d,q,xmean,var,trvar);
    %M-step
    k=feval('mstep',k,B,V,N,d,q,xmean,var,trvar);
    L1=feval('esteplike',k,B,V,N,d,q,xmean,var,trvar);
    if L1>L
      beep
      error('Che cazzo fai?\n');
         
    else
    end
    x=reshape(k(1:d*q),d,q);
    sigma=k(d*q+1);
    mu=k(d*q+2:end);
    for i=1:N %computes the inverses of the matrices M_n
      
      M(:,:,i)=inv(x'*F(:,:,i)*x+eye(q));
      xmean(:,i)=M(:,:,i)*x'*F(:,:,i)*(V(:,i)-mu');
      var(:,:,i)=M(:,:,i)+xmean(:,:,i)*xmean(:,:,i)';
      trvar(i)=trace(var(:,:,i));
      
      
    end 
   
  end
  w=k(1:d*q);
  sigma=k(d*q+1);
  mu=k(d*q+2:end);
  wML=reshape(w,d,q);
toc
beep
save results wML sigma mu
fprintf('Number of iterations used by EM algorythm\n')
t
fprintf('Output saved in results.mat\n')
