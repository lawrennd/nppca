function CNSPlotter_top50_3D(model)
%CNSPLOTTER_top50 plots results of PCA on CNS dataset on top 50 genes.



figure, plot3(model.W(1:7,1),model.W(1:7,2),model.W(1:7,3),'r+')
hold on
plot3(model.W(9:10,1),model.W(9:10,2),model.W(9:10,3),'r+')
hold on
plot3(model.W(34,1),model.W(34,2),model.W(34,3),'r+')
hold on
plot3(model.W(11:20,1),model.W(11:20,2),model.W(11:20,3),'bd')
hold on
plot3(model.W(31:33,1),model.W(31:33,2),model.W(31:33,3),'gx')
hold on
plot3(model.W(8,1),model.W(8,2),model.W(8,3),'gx')
hold on
plot3(model.W(35:end,1),model.W(35:end,2),model.W(35:end,3),'k>')
hold on
plot3(model.W(21:30,1),model.W(21:30,2),model.W(21:30,3),'mo')