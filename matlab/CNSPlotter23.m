function CNSPlotter23(model)
%CNSPLOTTER plots results of noisy PCA on CNS dataset

figure, plot(model.W(1:10,2),model.W(1:10,3),'r+')
hold on
 plot(model.W(11:20,2),model.W(11:20,3),'bd')
hold on
plot(model.W(21:24,2),model.W(21:24,3),'gx')
hold on
plot(model.W(25:32,2),model.W(25:32,3),'k>')
hold on
plot(model.W(33:end,2),model.W(33:end,3),'mo')