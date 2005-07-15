function [probes,Y, varY] = demNppcaLoadData(dataset)

% NPPCALOADDATA Load a dataset for demonstrating noisy PPCA.

% NPPCA


switch dataset
 
 case 'OC1_20'
   [probes, alpha, annotation, geneName] = gmosReadTxt(['../../gMOS/' ...
                       'data/signalOC1A.txt']);
   
   alpha = alpha(1:20, :);
   alphaLim = 1e-6;
   
   alpha(find(alpha<1e-6)) = 1e-6;
   
   Y = digamma(alpha);
   varY = trigamma(alpha);

 case 'OC1'
   [probes, alpha, annotation, geneName] = gmosReadTxt(['../../gMOS/' ...
                       'data/signalOC1B.txt']);
   
   alphaLim = 1e-6;
   alpha = alpha(5501:6000,:);
   probes= probes(5501:6000);

   alpha(find(alpha<1e-6)) = 1e-6;
   
   Y = digamma(alpha);
   varY = trigamma(alpha);

 case 'OC1B'
   [probes, alpha, annotation, geneName] = gmosReadTxt(['../../gMOS/' ...
                       'data/signalOC1B.txt']);
   
   alphaLim = 1e-6;
   alpha = alpha(5501:6000,:);
   probes= probes(5501:6000);

   alpha(find(alpha<1e-6)) = 1e-6;
   
   Y = digamma(alpha);
   varY = trigamma(alpha);

end