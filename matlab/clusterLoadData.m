function [annotation, Y, varY] = clusterLoadData()

% NPPCALOADDATA Load a dataset for demonstrating noisy PPCA.

% NPPCA
 

   [Y, annotation, geneName] = ClustergmosReadTxt(['../../gMOS/data/OC1_cluster_profile_nrec2.txt']);
   [varY, annotation, geneName] = ClustergmosReadTxt(['../../gMOS/data/OC1_cluster_Var2.txt']);
   
   Y = Y(1:2:199,:);
   varY = varY(1:2:199,:);
   annotation = annotation(1:2:199,:);

 