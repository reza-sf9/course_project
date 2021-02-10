clc; 
close all;
clear;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Prob1
syms w
%%% mexican hat
a=2;
f=sqrt(a)*((a*w)^2*sqrt(2*pi)*exp(-.5*((a*w)^2)));
ezplot(f,[-10 10 0 3]);
titl1=sprintf('frequency response of mexican hat wavelet with sigma=1 and a=2');
title(titl1)
grid on;


% %%% haar % 
% a=1;
% f=(sqrt(a))*(.5*(sin((a*w)/4)/((a*w)/4)));
% ezplot(f,[-30*pi 30*pi]);
% title('frequency response of Haar Wavelet, a=2 ');
% grid on;
n



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%  Prob 2

%%%%%%%%%%%%%%%%%%%%%%%%%%  caculate mexihat conv witha spcified a

b_step=.01; b_d=0; b_u=18;
b=b_d:b_step:b_u-b_step;


% f_b=1*(b>=6 & b<=12);
% 
% sigma=1; sigma_2=sigma^2;
% a=.25;
% 
% mex_d=-3.03*a; mex_u=3.03*a;
% b_m=mex_d:b_step:mex_u;
% mex1=(-1/(sigma_2*sqrt(a)))*(1-((b_m/a).^2)/sigma_2).*exp(-(((b_m/a).^2)/(2*sigma_2)));
% 
% f=conv(f_b,mex1);
% t=mex_d+b_d:b_step:b_u+mex_u-b_step;
% figure;
% plot(t,f)
% grid on;


%%%%%%%%%%%%%    calculate a  step by step with step a mexihat  
b_step=.01; b_d=0; b_u=18;
b=b_d:b_step:b_u-b_step;
f_b=1*(b>=6 & b<=12);

sigma=1; sigma_2=sigma^2;

mex_d=-3.03*1; mex_u=3.03*1;
b_m=mex_d:b_step:mex_u;
p=1;
for a=.25:.01:1

mex1=(-1/(sigma_2*sqrt(a)))*(1-((b_m/a).^2)/sigma_2).*exp(-(((b_m/a).^2)/(2*sigma_2)));
f(p,:)=conv(f_b,mex1);


p=p+1;
end
t=mex_d+b_d:b_step:b_u+mex_u-b_step;

figure;
plot(t,f(1,:));

aa=.25:.01:1;
[B,A]=meshgrid(t,aa);
mesh(B,A,f);
xlabel('b');ylabel('a');zlabel('WT(a,b)')



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Prob 3

syms w;
%%%%% Haaaar
phi_H=symfun(1/sqrt(2),w);
for k=1:1
   h=symfun(1/sqrt(2)*(1+exp(-1j*w/(2^k))),w);
   phi_H=1/sqrt(2)*phi_H*h;
end
figure
ezplot(abs(phi_H),[-1*pi 1*pi])
grid on;



Saw_H=symfun(1/sqrt(2),w);
g= 1/sqrt(2)*(1-exp(-1j*w/(2)));
Saw_H=Saw_H*g;
for k=2:100
   h=symfun(1/sqrt(2)*(1+exp(-1j*w/(2^k))),w);
   Saw_H=1/sqrt(2)*Saw_H*h;
end
figure;
ezplot(abs(Saw_H),[-15*pi 15*pi])
grid on;


%%% db2
phi_db=symfun(1/sqrt(2),w);
for k=1:10
   h=symfun(1/(4*sqrt(2))*((1+sqrt(3)) + (3+sqrt(3))*exp(-1j*w/(2^k)) + (3-sqrt(3))*exp(-1j*2*w/(2^k)) + (1-sqrt(3))*exp(-1j*3*w/(2^k))),w);
   phi_db=1/sqrt(2)*phi_db*h;
end
figure
ezplot(abs(phi_db),[-15*pi 15*pi 0 1])
grid on;



Saw_db=symfun(1/sqrt(2),w);
g= 1/(4*sqrt(2))*(1-exp(-1j*w/(2)));
Saw_db=Saw_db*g;
for k=2:100
   h=symfun(1/sqrt(2)*(1+exp(-1j*w/(2^k))),w);
   Saw_db=1/sqrt(2)*Saw_db*h;
end
figure;
ezplot(abs(Saw_db),[-15*pi 15*pi])
grid on;






% f=1;
% h=(1/sqrt(2))*[1 1];
% u = conv(f,h)
% 
% f=upsample(u,2)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% %%% Prob 4
% 
% syms f;
% g=symfun(0 ,f);
% for i=-100:100;
% %     g=g+(abs(sin(pi*(f+i)))^2)/abs(pi*(f+i))^2;
%       g=g+sinc(f+i)^2;
% end
% 
% ezplot(g,[-2 2 0 1.5])






%%%%%%%%%%%%%%%%%%%%%  Mexihat
% lb = -5;
% ub = 5;
% N = 100;
% [psi,xval] = mexihat(lb,ub,N);
% psi=-psi;
% figure;
% subplot(311);
% plot(xval,psi);
% title('Mexican Hat Wavelet');