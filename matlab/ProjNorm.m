function f=ProjNorm(X,U,D);
%PROJNORM computes the norm of a structured KL minimiser
d=diag(D);
uKL=X(1:end-1);
h=X(end);
f=-d'*((U'*[uKL';zeros(50,1)]).^2)-h*(uKL*uKL'-1);