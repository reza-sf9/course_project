clc; 
close all;
clear;


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%% Prob1
syms w;

%%% mexican hat
% a=2;
% f=sqrt(a)*((a*w)^2*sqrt(2*pi)*exp(-.5*((a*w)^2)));
% figure;
% ezplot(f,[-10 10 0 3]);
% titl1=sprintf('frequency response of mexican hat wavelet with sigma=1 and a=2');
% title(titl1)
% grid on;


%%% haar % 
a=1;
f=(sqrt(a))*(.5*(sin((a*w)/4)/((a*w)/4)));
figure;
ezplot(f,[-15*pi 15*pi]);
title('frequency response of Haar Wavelet, a=1 ');
grid on;

% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%  Prob 2
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%  caculate mexihat conv witha spcified a
% 
% b_step=.01; b_d=0; b_u=18;
% b=b_d:b_step:b_u-b_step;
% f=[];
% 
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
% 
% 
%%%%%%%%%%%%%    calculate a  step by step with step a mexihat  
b_step=.01; b_d=0; b_u=18;
b=b_d:b_step:b_u-b_step;
f_b=1*(b>=6 & b<=12);
f=[];

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

% 
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%% Prob 3
% f=[];
% 
% %%% part a & b
% %%%% choose your wavelet type and choose its function
% 
% %  w_type='haar';
%   w_type='db2';
% 
%   f_type='scaling';   %%  phi
% %   f_type='wavelet';   %%  saw
% %%%%%%
% 
% f=1;
% switch w_type
%     
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  Haar    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     case 'haar'
%         h=(1/sqrt(2))*[1 1];
%         g=(1/sqrt(2))*[1 -1];
%         
%         switch f_type
%             case 'scaling'   %% phi
%                 
%                 p=0;
%                 while p<10
%                     u = conv(f,h);
%                     f=upsample(u,2);
%                     p=p+1;
%                 end
%                 step=2/length(f);
%                 t=0:step:2-step;
%                 plot(t,f);
%                 title('scaling function (\phi) for Haar Wavelet')
%                 
%             case 'wavelet'  %% saw
%                 
%                 u=(1/sqrt(2))*conv(f,g);
%                 f=upsample(u,2);
%                 
%                 p=0;
%                 while p<9
%                     u = conv(f,h);
%                     f=upsample(u,2);
%                     p=p+1;
%                 end
%                 step=2/length(f);
%                 t=0:step:2-step;
%                 plot(t,f);
%                 title('scaling function (\psi) for Haar Wavelet')
%                 
%         end
%         
%         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  db2    %%%%%%%%%%%%%%%%%%%%%%%%%
%     case 'db2'
%         h=1/(4*sqrt(2))*[ 1+sqrt(3)    3+sqrt(3)  3-sqrt(3)   1-sqrt(3)];
%         g=1/(4*sqrt(2))*[ 1-sqrt(3)  -(3-sqrt(3)) 3+sqrt(3) -(1+sqrt(3))] ;
%         
%         switch f_type
%             case 'scaling'  %% phi
%                 
%                 p=0;
%                 while p<10
%                     u = conv(f,h);
%                     f=upsample(u,2);
%                     p=p+1;
%                 end
%                 step=2/length(f);
%                 t=0:step:2-step;
%                 plot(t,f);
%                 title('scaling function (\phi) for db2 Wavelet')
%                 
%                 
%             case 'wavelet'  %% saw
%                 
%                 u=(1/sqrt(2))*conv(f,g);
%                 f=upsample(u,2);
%                 
%                 p=0;
%                 while p<9
%                     u = conv(f,h);
%                     f=upsample(u,2);
%                     p=p+1;
%                 end
%                 step=2/length(f);
%                 t=0:step:2-step;
%                 plot(t,f);
%                 title('scaling function (\psi) for db2 Wavelet')
%         end
% end


%%% part c
w=[];
syms w;

%%%%%%% Haar
H_Haar = symfun(1/sqrt(2)*(1+exp(-1j*w)),w);

G_Haar = symfun(1/sqrt(2)*(1-exp(-1j*w)),w);




%%%%%%%% db2
H_db2 = 1/(4*sqrt(2))*((1+sqrt(3)) + (3+sqrt(3))*exp(-1j*w) + (3-sqrt(3))*exp(-1j*2*w) + (1-sqrt(3))*exp(-1j*3*w));
G_db2 = 1/(4*sqrt(2))*((1-sqrt(3)) - (3-sqrt(3))*exp(-1j*w) + (3+sqrt(3))*exp(-1j*2*w) - (1+sqrt(3))*exp(-1j*3*w));



figure;
ez1_h=ezplot(w,H_Haar,[0 pi]); grid on;
hold on
ez1_db2=ezplot(w,H_db2,[0 pi]); 
legend('H(\omega) (Haar)','H(\omega) (db2)'); 
title('Compare H(\omega) of Haar and H(\omega) of db2');


figure;
ez2_h=ezplot(w,G_Haar,[0 pi]); grid on;
grid on;
hold on
ez2_db2=ezplot(w,G_db2,[0 pi]); grid on;
legend('G(\omega) (Haar)','G(\omega) (db2)'); 
title('Compare G(\omega) of Haar and G(\omega) of db2');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Prob 4
f=[];
syms f;
g=symfun(0 ,f);
for i=-100:100;
%     g=g+(abs(sin(pi*(f+i)))^2)/abs(pi*(f+i))^2;
      g=g+sinc(f+i)^2;
end

ezplot(g,[-2 2 0 1.5])
