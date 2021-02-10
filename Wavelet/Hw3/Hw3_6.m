clc;
clear;
close all;

%%%%%%%%%%%%%%%%%%%%
Ts=1/400;                     %% sample rate
t=0:Ts:1;

x1=cos(2*pi*40*t);
x2=cos(2*pi*80*t);

% x2=x2+.5*randn(1,length(x2));

x=[x1 x2];
% x=exp(1i*.5*100*pi*t);

tprim=t.^2;
L=80;
W_L=2*L+1;                      %%window Length
DFT_L=8*W_L;                    %%DFT length
SG_L=length(x);

win=hanning(W_L);
x=[zeros(1,L),x,zeros(1,L)];
%%SG_L=length(x);

w=exp(1i*4*pi/DFT_L);
w=w.^(-L:L);                    %% length of exp must be equal to widow length


%%%%%%%%%%%%%%%%%%%%%%%    computing WVT    %%%%%%%%%%%%%%%%%%%%%%%%

B_F_WVT=DFT_L/4;                    %% we want to map until pi
Y_WVT=zeros(B_F_WVT,SG_L-1);
for n=1:SG_L-1
    s=x(n:n+W_L-1);
    s=s.*win';
    sf=fliplr(s);
    s=(sf.')'.*s;
    
    for k=0:B_F_WVT-1               
        r=(w.^k)';
        Y_WVT(k+1,n)=s*r;
    end
    
end
abs_Y_WVT=abs(Y_WVT);

%plot
Omega_WVT=(4*pi/DFT_L).*(1:B_F_WVT);
Time_WVT=(1:SG_L-1)/SG_L;
[time_WVT,omega_WVT]=meshgrid(Time_WVT,Omega_WVT);

figure;
mesh(time_WVT,omega_WVT,abs_Y_WVT);
xlabel('time');ylabel('frequency');zlabel('|WVT(n,k)|')

%%%%%%%%%%%%%%%%%%%%%%%    computing DSTFT    %%%%%%%%%%%%%%%%%%%%%%%%
% 
% w=exp(1i*2*pi/DFT_L);
% B_F_DSTFT=2*DFT_L;                          %% we want to map until pi
% Y_DSTFT=zeros(B_F_DSTFT,SG_L-W_L-1);      %% size of DSTFT
% 
% for n=1:SG_L-W_L-1
%     xx=x(n:n+W_L-1) .* win.';              %% calculates x[m]*w[n-m]
%     W=w.^(n-1:n+W_L-2);                     %% caculates exp
%     
%     for k=0:B_F_DSTFT-1                   %% this for samples from ASTFT
%         Y_DSTFT(k+1,n)=xx*(W.^k)';
%     end
% end
% 
% abs_Y_DSTFT=abs(Y_DSTFT);
% 
% Omega_DSTFT=2*pi*(0:B_F_DSTFT-1)/DFT_L;
% Time_DSTFT=1:SG_L-W_L-1;
% [time_DSTFT,omega_DSTFT]=meshgrid(Time_DSTFT,Omega_DSTFT);
% 
% 
% %%%%%%%%%%%%%%%%%%%%%% ploting WVT & DSTFT togther    %%%%%%%%%%%%%%%%%%%%
% figure;
% %%% plot WVT
% subplot(1,2,1);
% mesh(time_WVT,omega_WVT,abs_Y_WVT)
% xlabel('Time'),ylabel('Frequency'),zlabel('|WVT(n,k)|')
% title('WVT');
% 
% %%% plot DSTFT
% subplot(1,2,2);
% mesh(time_DSTFT,omega_DSTFT,abs_Y_DSTFT)
% xlabel('Time'),ylabel('Frequency'),zlabel('|DSTFT(n,k)|')
% title('DSTFT');