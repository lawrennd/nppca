

%gMOSPPCA Given the output of a gMOS reading, 
%performs a PPCA-like maximum likelihood optimisation on a latent variable 
%model.It expects to read from a file a vector of means, V, and a vector of 
%variances, beta (both row vectors).

  R=input('Which file should I read?','s');
     load(R);
     %Reads the file where the (column) vectors of means and variances are stored.
g=size(V);
N=g(2);  %N is the number of readings, usually ~10.
  d=g(1);  %d is the dimensionality of the readings, usually ~10000
h=min(B(:));
B=B./h;
V=V./sqrt(h);
    model.mu=mean(V'); %initialises the mean vector to the empirical mean
%Initialises the parameters of the LVM.
q=input('How many directions do you want to retain?\n');

tic;
     model.sigma=0.1;

     options = foptions;
%options(9) = 1;
options(1) = 1;
options(2)=1e-16;
options(3)=1e-16;
options(14)=5000;
%Initialises and starts the optimisation
G=cov((V-model.mu'*ones(1,N))');

[s,z]=pca(G,q);
%I use pca to initialise the factors;
model.W=z(:)';     
k=[model.W,model.sigma,model.mu]; %defines long vector with all parameters.
k= scg('gMOSPPCAlike',k,options,'gMOSPPCAgradlike',N,V,d,q,B);

model.W=k(1:d*q);
model.sigma=k(d*q+1)*sqrt(h);
model.mu=k(d*q+2:end)*sqrt(h);
wML=reshape(model.W,d,q)*sqrt(h);
toc
beep
save results wML model.sigma model.mu
fprintf('Output saved in results.mat\n')
%returns optimised parameters
