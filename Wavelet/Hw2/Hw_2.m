clc
close all;
clear;


t=0:1/150:1;     %%our nywueist rate is 130 
k=1;
for f=20:20:60
    x(:,k)=(cos(2*pi*f*t)+cos(2*pi*(f+5)*t));
    k=k+1;
end

x=x(:);%%put 3 rows in 1 row
x=x.';
Lx=length(x);
noise=.5*randn(1,Lx);
%x=x+noise;
L=60;                    %window length
N=length(x);

% x1=x(1:L);
% x1=fliplr(x1);
% x2=x(N-L:N);
% x2=fliplr(x2);


%zero pad k aval va entehaa panjere daghighan rooye signal gharar girad
x=[zeros(1,L),x,zeros(1,L)];
figure;
plot(x);
N=length(x);
NN=6*L;                     %length of dft
w=exp(1i*2*pi/NN);
window=hanning(L);
window=window.';

piS=NN;                        % if piS = DFT_L , we samples until 2*pi 
                               % if we want sample less than 2*pi we must set the
                               % current value for piS
%stft
for m=1:N-L
    xx=x(m:m+L-1).*window;     % calculates x[m]*w[n-m]
    W=w.^(m-1:m+L-2);          % caculates exp
    for k=0:piS-1
        y(m,k+1)=xx*(W.^k)';   % ' operator perform transpose and conjucate 
    end
end
yy=abs(y);
omega=2*pi*(0:piS-1)/NN;
T=0:N-L-1;
[omega,time]=meshgrid(omega,T);
mesh(time,omega,yy);





