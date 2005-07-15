%DEMKLPCA2 demo showing structured KL minimisation vs block diagonal PCA.
dim = 2*50;
randn('seed', dim);
a=randn(dim);
G=a*a';
G1 = G(1:dim/2,1:dim/2);
[var1,u1,lambda1]= ppca(G1,1);
[U,D]=eig(G);
options=foptions;
options(9)=1;
options(10)=1;
%options(14)=1000;
uKL=[1,zeros(1,49)];
h=1;
X=[uKL,h]; %Introduce vector containing best eigenvalue and Lagrange multiplier
X=scg('ProjNorm', X, options, 'GradProjNorm',U,D);
uKL=X(1:end-1);
plot(uKL,u1,'r+');
%axis equal