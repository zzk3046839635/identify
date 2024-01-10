% LQ decomposition
function [L11,L21,L22,L31,L32,L33,Q1T,Q2T,Q3T]=lqBCA(B,C,A)
kB=size(B,1);
kC=size(C,1);
[Q,R]=qr([B;C;A]');
QT=Q';
L=R';%[B;C;A]=L*QT
L11=L(1:kB,1:kB);
L21=L(kB+1:kB+kC,1:kB);
L22=L(kB+1:kB+kC,kB+1:kB+kC);
sLr=size(L,1);
sLc=size(L,2);
L31=L(kB+kC+1:sLr,1:kB);
L32=L(kB+kC+1:sLr,kB+1:kB+kC);
L33=L(kB+kC+1:sLr,kB+kC+1:sLc);
Q1T=QT(1:kB,:);
Q2T=QT(kB+1:kB+kC,:);
sQTr=size(QT,1);
Q3T=QT(kB+kC+1:sQTr,:);

