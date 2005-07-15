%DEMKLPCA demo showing structured KL minimisation vs block diagonal PCA.
a=[1,1,0.3,0.3];
b=[0.3,0.3,2,1];
G=a'*a+b'*b+0.1*eye(4);
G1 = G(1:2,1:2);
G2 = G(3:4,3:4);
[var,U,lambda] = ppca(G,2);
[var1,u1,lambda1]= ppca(G1,1);
[var2,u2,lambda2]= ppca(G2,1);
U2=[u1;0;0];
U1=[0;0;u2];
UKL1=[0;0;U(3:4,1)]/sqrt(U(3,1)^2+U(4,1)^2);
UKL2=[U(1:2,2);0;0]/sqrt(U(1,2)^2+U(2,2)^2);
UKL=[UKL1,UKL2];
x=gsamp(zeros(1,4),G,100);
Proj=x*U;
plot(Proj(:,1),Proj(:,2),'r+');
hold on
axis equal
ProjBlock=[x*U1,x*U2];
plot(ProjBlock(:,1),ProjBlock(:,2),'bd');
hold on
ProjKL=x*UKL;
plot(ProjKL(:,1),ProjKL(:,2),'gx');
