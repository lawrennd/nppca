function []=demProfilePlotter(data,vcf)

%DEMPROFILEPLOTTER demo showing comparison of reconstructed and
%original profile for two genes.

set(gca, 'fontsize', 20)
set(gca, 'fontname', 'helvetica')

days = [0:9 11 14];
switch data
  
 case 'gata3'
  
  if vcf==1
    load('risultatigata3')
  elseif vcf==5
      load risultatigata3vcf5
  else
    load risultatigata3vcf10
  end
  
  [probes,annotation,Y,varY] = nppcaLoadData('OC1B');
  probes=probes(5501:6000);
  annotation=annotation(5501:6000);
  Y=Y(5501:6000,:);
  varY=varY(5501:6000,:)/vcf;
  figure, z = errorbar(days,Y(454,:), sqrt(varY(454,:)));
  
  [S, varS]=reconstruct(model, expectations);
  hold on
  
  h = errorbar(days,S(454,:), sqrt(varS(454,:)),'r-');
  set(h, 'linewidth', 2)
  set(h, 'linestyle', '--')
  hold off
  ylabel('expression levels');
 xlabel('days');
 
 otherwise
  
  if vcf==1
     load risultatip27
   elseif vcf==5
      load risultatip27vcf5
   else
    load risultatip27vcf10
   end
   
   [probes,annotation,Y,varY] = nppcaLoadData('OC1B');
   probes=probes(3501:4000);
annotation=annotation(3501:4000);
Y=Y(3501:4000,:);
varY=varY(3501:4000,:)/vcf;
figure, z = errorbar(days,Y(160,:), sqrt(varY(160,:)));

[S, varS]=reconstruct(model, expectations);
  
  figure, h = errorbar(days,S(160,:), sqrt(varS(160,:)),'r-');
  set(h, 'linewidth', 2)
  set(h, 'linestyle', '--')
  hold off
  ylabel('expression levels');
 xlabel('days');
end
