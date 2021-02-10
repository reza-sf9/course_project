clc;
clear;
close all;

eta1=6; t=.69;
a=(eta1-t)/sqrt(2*eta1);
b=(eta1+t)/sqrt(2*eta1);

c=.6*(qfunc(a)) + .4*(qfunc(b))


kj


%%%%%%%%%%%%%%%%%%%%
Ts=1/400;                          %% sample rate
t=0:Ts:1;

x1=cos(2*pi*20*t);
x2=cos(2*pi*40*t);

x2=x2+.5*randn(1,length(x2));

x=[x1 x2];
figure;
plot(x);
title('original signal');


W_L=80;                            %% Window Length
window=hanning(W_L);   
window=window.';

%%% at the first we must do zero-pad ,cause when we use window ,
%%% window should put on the first sample of signal
x=[zeros(1,W_L),x,zeros(1,W_L)];
SG_L=length(x);                     %% Signal Length
DFT_L=4*W_L;                        %% DFT Length
w=exp(1i*2*pi/DFT_L);


%%% DSTFT
B_F=DFT_L;                          %% (Band Frequency)this variable detemines bound of frequeny ,that shown by DSTFT 
                                    %%% (indeed deteminies ,where is the bound of our sampling from ASTFT )
                                    
y=zeros(B_F,SG_L-W_L-1);            %% size of DSTFT

for m=1:SG_L-W_L-1
    xx=x(m:m+W_L-1).*window;        %% calculates x[m]*w[n-m]
    W=w.^(m-1:m+W_L-2);             %% caculates exp
    
    for k=0:B_F-1                   %% this for samples from ASTFT
        y(k+1,m)=xx*(W.^k)';        
    end
end

thr=8.3;
y_S_Thr=wthresh(y,'s',thr);            %% soft thresholding
% Abs_y_Sthr=abs(y_Sthr)+thr;           %% this term comensate the reduced value from X(n,k) by soft thresholding 
% Angle_y_Sthr=angle(y_Sthr);
% y_n_Sthr=Abs_y_Sthr.*(exp(1i).^Angle_y_Sthr);


Orig_S=abs(y);                       %% DSTFT of Original Signal
DeN_S=abs(y_S_Thr);               %% DSTFT of Denoised Signal

Omega=2*pi*(0:B_F-1)/DFT_L;
T=1:SG_L-W_L-1;
[time,omega]=meshgrid(T,Omega);
figure;

%%% plot DSTFT of Original Signal

subplot(1,2,1);
mesh(time,omega,Orig_S)
xlabel('Time'),ylabel('Frequency'),zlabel('|F(n,k)|')
title('DSTFT of Signal with Noise');

%%% plot DSTFT of Denoised Signal
subplot(1,2,2);
mesh(time,omega,DeN_S)
xlabel('Time'),ylabel('Frequency'),zlabel('|F(n,k)|')
title('DSTFT of Signal after Soft Thresholding');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%% Implemntation of Inverse DSTFT %%%%%%%%%%%%%

%%%%%%%%%%%   filter bank sum (FBS) method   %%%%%%%%%%%%%         %% w=80 
sig=y_S_Thr; 

[k,n]=size(sig);
N=DFT_L;
x_FBS=zeros(n,1);
w0=window(1);

for r1=0:n-1  
    
    temp1=(exp(1i*2*pi*r1/N).^(0:k-1))*sig(:,r1+1);
    
    x_FBS(r1+1)=(1/(w0*N)).*temp1;
end

 x_FBS=real(x_FBS(W_L:end));
 %%% scaling
 a=abs(x_FBS(50:350));
 b=max(a);
 x_FBS=x_FBS/b;

figure;
plot(x_FBS);
title('ISTFT (FBS approach)');


figure;
subplot(1,2,1);
plot(x);
title('original signal');
subplot(1,2,2);
plot(real(x_FBS(W_L:end)));
title('ISTFT (FBS approach)');
%%%%%%%%%%%%%%%%    overlap add (OLA) method   %%%%%%%%%%%%%    w=250

x_OLA=zeros(n,1);
W_S = sum(window);       % Summation of Window
 
for nn=0:n-1
        temp3=0;
 for    r2=0:n-1
        temp2=(exp(1i*2*pi*nn/N).^(0:k-1))*sig(:,r2+1);
        
        temp3=temp3+(1/N).*temp2;
 end
   x_OLA(nn+1)=temp3;
end

a=abs(x_OLA(100:200));
b=max(a);
% x_OLA=x_OLA/b;


 x_OLA = real(x_OLA/W_S);
 
 %%% scaling
 a=abs(x_OLA(50:350));
 b=max(a);
 x_OLA=x_OLA/b;
 
 x_OLA=x_OLA(W_L:end);
 figure;
 plot(x_OLA);
 title('ISTFT (OLA approach)');
