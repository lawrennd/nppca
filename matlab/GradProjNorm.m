function f=GradProjNorm(X,U,D)
%GRADPROJNORM Gradient of ProjNorm
d=diag(D);
uKL=X(1:end-1);
h=X(end);
T=(d.*(U'*[uKL';zeros(50,1)]))'*U';
S=T(1:50);
f=-[2*S-2*h*uKL,(norm(uKL)-1)];