function CNSPlotter3D(model)
%CNSPLOTTER plots results of noisy PCA on CNS dataset

figure, plot3(model.W(1:10,1),model.W(1:10,2),model.W(1:10,3),'r+')
hold on
 plot3(model.W(11:20,1),model.W(11:20,2),model.W(11:20,3),'bd')
hold on
plot3(model.W(21:24,1),model.W(21:24,2),model.W(21:24,3),'gx')
hold on
plot3(model.W(25:32,1),model.W(25:32,2),model.W(25:32,3),'k>')
hold on
plot3(model.W(33:end,1),model.W(33:end,2),model.W(33:end,3),'mo')