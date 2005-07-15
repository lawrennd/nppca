%PROFILECOMPARE Compares profiles of GATA3, P27 and mean profile

load resultsGata3
modelGata=model;
expectationsGata=expectations;

load resultsp27
modelP27=model;
expectationsP27=expectations;

days = [0:9 11 14];

[SGata, varSGata]=reconstruct(modelGata, expectationsGata);

figure, h = errorbar(days,SGata(454,:), sqrt(varSGata(454,:)),'r*--');

hold on
set(h, 'linewidth', 2)


plot(days, modelGata.mu, 'bo-');
hold on



[Sp27,varSp27]=reconstruct(modelP27,expectationsP27);

z = errorbar(days,Sp27(160,:), sqrt(varSp27(160,:)),'md-.');

hold on
set(z, 'linewidth', 2)


plot(days,modelP27.mu,'c+:');
set(gca, 'fontsize', 20)
set(gca, 'fontname', 'helvetica')
set(gca,'ylim', [-2 5])
set(gca, 'xlim', [-0.3 14]);
set(gca, 'xtick', [0:2:14])
hold off
ylabel(['log expression levels']);
xlabel('days');