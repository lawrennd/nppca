function f=CNSPlotAssesser(model,nComp);

%PLOTASSESSER Assesses the clustering of PCAor NPPCA using a 1NN
%classifier for CNS dataset.
npts=size(model.W,1);
X=model.W(:,1:nComp);
Distances=dist2(X,X);
Distances=Distances-diag(diag(Distances))+diag(100*ones(1, npts));
[Min,I]=min(Distances);
score=zeros(1,npts);
for i=1:10
  if I(i)>10
    score(i)=1;
  else
  end
end
for i=11:20
  if I(i)>20
    score(i)=1;
  elseif I(i) < 11
    score(i)=1;
  else
  end
end
for i=21:24
  if I(i)>24
    score(i)=1;
  elseif I(i) < 21
    score(i)=1;
  else
  end
end
for i=25:32
  if I(i)>32
    score(i)=1;
  elseif I(i) < 25
    score(i)=1;
  else
  end
end
for i=33:42
  if I(i)<33
    score(i)=1;
  
  else
  end
end
f=sum(score)/npts;