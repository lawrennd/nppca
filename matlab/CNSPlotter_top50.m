function CNSPlotter_top50(model)
%CNSPLOTTER plots results of noisy PCA on CNS dataset

figure, plot(model.W(1:10,1),model.W(1:10,2),'r+')
hold on
 plot(model.W(11:20,1),model.W(11:20,2),'bd')
hold on
plot(model.W(31:34,1),model.W(31:34,2),'gx')
hold on
plot(model.W(35:end,1),model.W(35:end,2),'k>')
hold on
plot(model.W(21:30,1),model.W(21:30,2),'mo')