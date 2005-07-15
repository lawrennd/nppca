function nppcaProfilePlotter(model,expectations,Y,varY,ngene,geneName)

%NPPCAPROFILEPLOTTER shows comparison of reconstructed and
%original profile for genes.



days = [0:9 11 14];

figure, z = errorbar(days,Y(ngene,:), sqrt(varY(ngene,:)));

[S, varS]=reconstruct(model, expectations);
hold on

h = errorbar(days,S(ngene,:), sqrt(varS(ngene,:)),'r-');
set(h, 'linewidth', 2)
set(h, 'linestyle', '--')
set(gca, 'fontsize', 20)
set(gca, 'fontname', 'helvetica')
%set(gca,'ylim', [-5 5])
set(gca, 'xlim', [-0.3 14]);
set(gca, 'xtick', [0:2:14])
hold off
ylabel(['log expression levels']);
xlabel('days');


