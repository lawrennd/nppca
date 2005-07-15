%CLUSTEROC1 clustering for output of OC1 experiment
clear functions
clear all


load demOC1I;
[probesA, annotationsA,YA, varYA]=nppcaLoadData('OC1A');
[probesB, annotationsB,YB, varYB]=nppcaLoadData('OC1B');
Y=[YA;YB];
varY=[varYA;varYB];
probes=[probesA;probesB];
annotations=[annotationsA;annotationsB];

nPC=size(model.W,2);
expectations.x=zeros(size(Y,1),nPC);
expectations=nppcaEstep(model, expectations, varY, Y);

[exp_sorted,I]=sort(expectations.x);

L=100; %number of genes selected from the top 
[Recon_prof, Recon_var]=reconstruct(model, expectations);
for j=1:nPC
  cluster.an(:,j)=annotations(I(end:-1:end-L+1,j));
  cluster.probes(:,j)=probes(I(end:-1:end-L+1,j));
  cluster.profile(:,:,j)=Recon_prof(I(end:-1:end-L+1,j),:);
  cluster.var(:,:,j)=Recon_var(I(end:-1:end-L+1,j),:);


  outputFile=['OC1_cluster_profile',num2str(j),'.txt'];
  f = fopen(outputFile,'w');
  fprintf('Writing: %s',char(outputFile))
  fprintf(f,'%s\t %s\t %s\t %s\t %s\t %s\t %s\t %s\t %s\t %s\t %s\t %s\t %s\t %s\t \n',...
          'Probe Name','Gene Name',...
          'Signal 0d', 'Signal 1d','Signal 2d','Signal 3d','Signal 4d',...
          'Signal 5d','Signal 6d','Signal 7d','Signal 8d','Signal 9d',...
          'Signal 11d','Signal 12d');
  
  
  for i=1:L
    fprintf(f,'%s\t %s\t %6.4f\t %6.4f\t %6.4f\t %6.4f\t %6.4f\t %6.4f\t %6.4f\t %6.4f\t %6.4f\t %6.4f\t %6.4f\t %6.4f\t \n',...
            char(cluster.probes(i,j)),char(cluster.an(i,j)),...
            cluster.profile(i,1,j),cluster.profile(i,2,j),cluster.profile(i,3,j),...
            cluster.profile(i,4,j),cluster.profile(i,5,j),cluster.profile(i,6,j),...
            cluster.profile(i,7,j),cluster.profile(i,8,j),cluster.profile(i,9,j),...
            cluster.profile(i,10,j),cluster.profile(i,11,j),cluster.profile(i,12,j));
  end
  fclose(f);
  fprintf(' Done \n\n')
  
  outputFile=['OC1_cluster_Var',num2str(j),'.txt'];
  f = fopen(outputFile,'w');
  fprintf('Writing: %s',char(outputFile))
  fprintf(f,'%s\t %s\t %s\t %s\t %s\t %s\t %s\t %s\t %s\t %s\t %s\t %s\t %s\t %s\t \n',...
          'Probe Name','Gene Name',...
          'Variance 0d', 'Variance 1d','Variance 2d','Variance 3d','Variance 4d',...
          'Variance 5d','Variance 6d','Variance 7d','Variance 8d','Variance 9d',...
          'Variance 11d','Variance 12d');
  
  
  for i=1:L
    fprintf(f,'%s\t %s\t %6.4f\t %6.4f\t %6.4f\t %6.4f\t %6.4f\t %6.4f\t %6.4f\t %6.4f\t %6.4f\t %6.4f\t %6.4f\t %6.4f\t \n',...
            char(cluster.probes(i,j)),char(cluster.an(i,j)),...
            cluster.var(i,1,j),cluster.var(i,2,j),cluster.var(i,3,j),...
            cluster.var(i,4,j),cluster.var(i,5,j),cluster.var(i,6,j),...
            cluster.var(i,7,j),cluster.var(i,8,j),cluster.var(i,9,j),...
            cluster.var(i,10,j),cluster.var(i,11,j),cluster.var(i,12,j));
  end
  fclose(f);
  fprintf(' Done \n\n')

end
