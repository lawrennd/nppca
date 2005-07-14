function f=CNSPlotAssesser_top50(model,nComp);

%PLOTASSESSER Assesses the clustering power of a plot using a 1NN
%classifier.
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
for i=21:30
  if I(i)>30
    score(i)=1;
  elseif I(i) < 21
    score(i)=1;
  else
  end
end
for i=31:34
  if I(i)>34
    score(i)=1;
  elseif I(i) < 31
    score(i)=1;
  else
  end
end
for i=35:42
  if I(i)<35
    score(i)=1;
  
  else
  end
end
f=sum(score)/npts;