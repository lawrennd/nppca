function [Probes, Annotation, Y, varY] = RALoadData(dataset)
% NPPCALOADDATA Load a dataset for demonstrating noisy PPCA.

% NPPCA

switch dataset
 
 case 'TOYRA_20'

   [Probes, Y, Annotation] = read_RAtxtFiles(['../../gMOS/data/RA_LogSignal.txt']);
   [Prob, varY, Annotat] = read_RAtxtFiles(['../../gMOS/data/RA_LogSignalVar.txt']);   
  Y = Y(1:20, :);
  Probes = Probes(1:20,:);
  Annotation = Annotation(1:20,:);
  varY=varY(1:20,:);
   YLim = -1e6;
   
   Y(find(Y<YLim)) = YLim;
    
    case 'RA'
        
        [Probes, Y, Annotation] = read_RAtxtFiles(['../../gMOS/' ...
                       'data/RA_LogSignal.txt']);
   [Probes, varY, Annotation] = read_RAtxtFiles(['../../gMOS/' ...
                       'data/RA_LogSignalVar.txt']);   

   YLim = -1e6;
   
   Y(find(Y<YLim)) = YLim;
end
