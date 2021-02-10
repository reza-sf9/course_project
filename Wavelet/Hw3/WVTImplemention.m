clc
clear
close all;

t=0:1/400:1;

%f=cos(2*pi*(5t+2+5*(t.^2)))
x1=cos(2*pi*50*t);
x2=cos(2*pi*75*t);

f=[x1 x2];
sig_length=length(f);
noise=.1*randn(1,sig_length);
f=f+noise;

L=40;
W_L=2*L+1;                      %%window Length
DFT_L=8*W_L;                    %%DFT length

win=boxcar(W_L);
f=[zeros(1,L),f,zeros(1,L)];
w=exp(1i*4*pi/DFT_L);
w=w.^(-L:L);                    %% length of exp is equal to widow length


for n=1:sig_length-1
    s=f(n:n+W_L-1);
    s=s.*win';
    sf=fliplr(s);
    s=(sf.')'.*s;

    for k=0:DFT_L/4%% we want to map until pi
        r=(w.^k)';
        F(n,k+1)=s*r;
    end
    
 end
    k=(4*pi/DFT_L).*(0:DFT_L/4);
    n=(1:sig_length-1);
    figure;
    [n,k]=meshgrid(n,k);
    mesh(n,k,abs(F).');
    xlabel('time');ylabel('frequency');zlabel('|F|')
    