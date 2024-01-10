clc
clear
load('input.mat')%input mat name
load('output.mat')%output mat name
l=1;
m=1;
NN=5000;%number of datas

i=95;%choose suitable value
j=3000;%choose suitable value

u=in.signals.values;
y=out1.signals.values(:,1)';

Yp=zeros(i*l,j);
Up=zeros(i*m,j);
Yf=zeros(i*l,j);
Uf=zeros(i*m,j);
for k1=1:i*l
    for k2=1:j
        Yp(k1,k2)=y(k1+k2-1);
        Yf(k1,k2)=y(i*l+k1+k2-1);
    end
end
for k1=1:i*l
    for k2=1:j
        Up(k1,k2)=u(k1+k2-1);
        Uf(k1,k2)=u(i*l+k1+k2-1);
    end
end
Wp=[Yp;Up];

[Oid,Oilq]=obp(Uf,Wp,Yf);
[UU,SS,VV]=svd(Oilq);%SVD
n=2;%system order
U1=UU(:,1:n);
S1=SS(1:n,1:n);
T=[1 2;0 2];%any matrix with system order
VVT=VV';
V1T=VVT(1:n,:);
Gammai=U1*sqrtm(S1)*T;
Xf=inv(T)*sqrtm(S1)*V1T;
Xip1jm1=Xf(:,2:j);
Xijm1=Xf(:,1:j-1);
Uijm1=zeros(m,j-1);
Yijm1=zeros(l,j-1);
for k1=1:j-1
    Uijm1(:,k1)=u(i+k1);
    Yijm1(:,k1)=y(i+k1);
end
hatABCD=[Xip1jm1;Yijm1]*[Xijm1;Uijm1]'*inv([Xijm1;Uijm1]*[Xijm1;Uijm1]');
AT=hatABCD(1:n,1:n)
BT=hatABCD(1:n,n+1:n+m) 
CT=hatABCD(n+1:n+l,1:n)
DT=hatABCD(n+1:n+l,n+1:n+m)

xT=zeros(2,NN+1);
for k=1:NN
    xT(:,k+1)=AT*xT(:,k)+BT*u(k);
    yT(k)=CT*xT(:,k)+DT*u(k);
    t(k)=k-1;
end
%show result
ltext={'y(k)','yT(k)'};
h(1,1)=plot(t(1:3000),y(1:3000),'-b','LineWidth',2);hold on;
h(2,1)=plot(t(1:3000),yT(1:3000),'--r','LineWidth',2);hold off;
legend(h,ltext);
xlabel('time sequence k')
ylabel('y(k) and yT(k)')

%[num2, den2]=ss2tf(AT, BT, CT, DT, 1);
%printsys(num2, den2);

%eig(AT)
ess=sqrt(sum((y(1:NN)-yT(1:NN)).^2)/sum(y(1:NN).^2)/NN)*100;
disp(ess);
