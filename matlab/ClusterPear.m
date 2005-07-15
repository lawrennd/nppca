%CLUSTERPEAR Performs a Pearson coefficient based hierarchical clustering
%for gMOS processed data, downweighting the high-variance genes.
[annotations, Y, varY]=clusterLoadData();
Ymean=mean(Y,1);

dataDim=size(Y,1);
numData=size(Y,2);
normalised_Y=zeros(dataDim,numData);
for i=1:dataDim
    normalised_Y(i,:)=(Y(i,:)-Ymean)/norm(Y(i,:)-Ymean);
    newVarY(i,:) = ones(1,numData)+varY(i,:)/norm(Y(i,:)-Ymean);
end
for k=1:dataDim
    Corr=zeros(dataDim-k+1,dataDim-k+1);
    for i=1:dataDim-k+1
        for j=i:dataDim-k+1
            Corr(i,j)=Prob(normalised_Y(i,:),normalised_Y(j,:),varY(i,:),varY(j,:));
            Corr(j,i)=Corr(i,j);
        end
        Corr(i,i)=0;
    end
    [max_col,I]=max(Corr);
    [B,J]=max(max_col); 
    Dist(k)=Corr(I(J),J);
    Cluster={annotations{I(J)};annotations{J}}';
    newMean=(normalised_Y(I(J),:).*(varY(I(J),:)).^(-1)+...
        (varY(J,:)).^(-1).*normalised_Y(J,:))./(varY(I(J),:).^(-1)+...
        varY(J,:).^-1);
    newNormalisedMean=newMean/norm(newMean);
    normalised_Y=[normalised_Y(1:min(I(J),J)-1,:);normalised_Y(min(I(J),J)+1:max(I(J),J)-1,:);normalised_Y(max(I(J),J)+1:end,:);newNormalisedMean];
    varY=[varY(1:min(I(J),J)-1,:);varY(min(I(J),J)+1:max(I(J),J)-1,:);varY(max(I(J),J)+1:end,:);2*(varY(I(J),:).*varY(J,:))./(varY(I(J),:)+varY(J,:))];
    annotations={annotations{1:min(I(J),J)-1},annotations{min(I(J),J)+1:max(I(J),J)-1},annotations{max(I(J),J)+1:end},Cluster}';
end
save resultCluster annotations