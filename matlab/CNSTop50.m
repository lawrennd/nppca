%CNSTOP50 recreates experiment of Pomeroy et al.
[genes, Y50]=readCNS_top50('../data/top_mas4.txt');
[PC50.var,PC50.W]=pca(Y50);
CNSPlotter_top50_3D(PC50);