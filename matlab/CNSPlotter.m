function CNSPlotter(model)
%CNSPLOTTER plots results of noisy PCA on CNS dataset

figure, plot(model.W(1:10,1),model.W(1:10,2),'r+')
hold on
 plot(model.W(11:20,1),model.W(11:20,2),'bd')
hold on
plot(model.W(21:24,1),model.W(21:24,2),'gx')
hold on
plot(model.W(25:32,1),model.W(25:32,2),'k>')
hold on
plot(model.W(33:end,1),model.W(33:end,2),'mo')