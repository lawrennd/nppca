function [S, varS]=reconstruct(model, expectations)

% RECONSTRUCT reproduces a denoised reconstruction of the gene expressions with errorbars.

% NPPCA

N=size(expectations.x,1);
d=size(model.W,1);
S=expectations.x*model.W'+ones(N,1)*model.mu;
varS=zeros(N, d);
for i=1:N
  varS(i,:)=diag(model.W*(expectations.xxT(:,:,i)- ...
                                 expectations.x(i,:)'*expectations.x(i,:))*model.W')';
  if varS(i)<0
    error('Negative variance')
  end
end
