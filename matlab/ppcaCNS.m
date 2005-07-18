% PPCA for the CNS data

[probes,Y, geneName] =readCNS_DATA('../data/Dataset_A_multiple_tumor_samples.txt');
[PcCoeff, PcVec] = pca(Y,3);
MD=PcVec(1:10,:);
MGlio=PcVec(11:20,:);
Rhab=PcVec(21:30,:);
NCer=PcVec(31:34,:);
PNET=PcVec(35:end,:);

figure, plot3(MGlio(:,1),MGlio(:,2),MGlio(:,3),'bd');
hold on
plot3(MD(:,1),MD(:,2),MD(:,3),'r+');
hold on
plot3(Rhab(:,1),Rhab(:,2),Rhab(:,3),'m>');
hold on
plot3(NCer(:,1),NCer(:,2),NCer(:,3),'go');
hold on
plot3(PNET(:,1),PNET(:,2),PNET(:,3),'c*');
grid on