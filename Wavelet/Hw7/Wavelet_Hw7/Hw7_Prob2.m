clc;
clear;
close all;


%%%%% Prob 2 %%%%%%%

step = 0.0005; f_d=0; f_u=2;
t_f = f_d : step : f_u ;
f=zeros(length(t_f),3);
for i=1:3
    landa=.1 * 2^i;
    f(:,i)= 1- abs(0.5 - t_f).^landa;
end
f1=f(:,1); f2=f(:,2); f3=f(:,3);

figure;
subplot(311), plot(t_f,f1), title('f1 ***  \lambda = 0.2');
subplot(312), plot(t_f,f2), title('f2 ***  \lambda = 0.4');
subplot(313), plot(t_f,f3), title('f3 ***  \lambda = 0.8');


Gu_d=-6; Gu_u=6;
t_Gu=Gu_d:step:Gu_u;         

t=Gu_d + f_d :step: Gu_u + f_u;

f1_conv=zeros(64,length(t));    f2_conv=zeros(64,length(t));    f3_conv=zeros(64,length(t));
f1_conv_d1=zeros(64,length(t)); f2_conv_d1=zeros(64,length(t)); f3_conv_d1=zeros(64,length(t));
max_WT_f1_d1=zeros(1,64);          max_WT_f2_d1=zeros(1,64);          max_WT_f3_d1=zeros(1,64);

for s=1:64

Gu=(sqrt(1)/sqrt(2*pi*s)).*exp(-(t_Gu/(sqrt(2))*s).^2);

%%%% f1
f1_conv(s,:) =    conv(f1,Gu);                                   %%  calculate f(u) * (sqrt()×Psi((-u*s))
f1_conv_d1(s,:) = (gradient(f1_conv(s,:))/step);                %% first derivative of f1


singular_index=find(f1_conv(s,:)==max(f1_conv(s,:)));        %% obtain the index's number of t= 0.5 
max_WT_f1_d1(s)= max(f1_conv_d1(s,(singular_index-500:singular_index+500)));  %% caclulate Local Max for t=.5                       %% caclulate Local Max for t=.5

%%%% f2
f2_conv(s,:) = conv(f2,Gu);
f2_conv_d1(s,:) = gradient(f2_conv(s,:))/step;          

singular_index=find(f2_conv(s,:)==max(f2_conv(s,:)));
max_WT_f2_d1(s)= max(f2_conv_d1(s,singular_index-500:singular_index+500));
% max_WT_f2_d1(s)= max(f2_conv_d1(s,:));

%%%% f3
f3_conv(s,:) = conv(f3,Gu);
f3_conv_d1(s,:) = gradient(f3_conv(s,:))/step;                

singular_index=find(f3_conv(s,:)==max(f3_conv(s,:)));
max_WT_f3_d1(s)= max(f3_conv_d1(s,singular_index-500:singular_index+500));
end 



%%%% Mesh plot

% for p=1:59
% temp(p,:)=f3_conv_d1(p+5,:);
% end
% figure;
% s=6:64;
% [U,S]=meshgrid(t,s);
% mesh(U,S,temp);
% xlabel('U');ylabel('S');zlabel('WT(U,S)  related to f3')


log_s = log(1:64);
log_WT_f1= -log(max_WT_f1_d1);
log_WT_f2= -log(max_WT_f2_d1);
log_WT_f3= -log(max_WT_f3_d1);

%%%% plotting Modulus Maxima Characteristi
figure;
plot((log_s(5:64)),log_WT_f2(5:64),'g*')
title('Modulus Maxima Characteristic for f3')


%%%% estimation lambda for f1
m_1=polyfit(log_s(10:64),log_WT_f1(10:64),1);
alpha_1=(m_1(1)-.5);
%%%% estimation lambda for f2
m_2=polyfit(log_s(10:64),log_WT_f2(10:64),1);
alpha_2=(m_2(1)-.5);
%%%% estimation lambda for f3
m_3=polyfit(log_s(10:64),log_WT_f3(10:64),1);
alpha_3=(m_3(1)-.5);
plk