clc
clear
close all

%%%%%%%%%%% prob 2
% c_0=1/sqrt(2).*[3 8 21 5 0 -18 -7 -3];
% d_0=1/sqrt(2).*[-1 -2 -5 1 2 -8 -1 -1];
% C=[c_0 d_0];
% L=[8 8 16];
% % a0 = waverec(c_0k,d_0k,'db1');
% 
% a1= waverec(C,L,'db1');

%%%%%%%%%%% prob 3
Ts=3/800;
t=(0:Ts:2);
f=cos(2*pi*50*t);  %% C_0,k

[c,l]=wavedec(f,4,'db1');

d{1}=c(l(1)+l(2)+l(3)+l(4)+1:l(1)+l(2)+l(3)+l(4)+l(5));
d{2}=c(l(1)+l(2)+l(3)+1:l(1)+l(2)+l(3)+l(4));
d{3}=c(l(1)+l(2)+1:l(1)+l(2)+l(3));
d{4}=c(l(1)+1:l(1)+l(2));
c_4=c(1:l(1));


% %%%% Plot Decomposition's Coeficents
figure;
subplot(5,1,1); plot(f); title('original signal')
subplot(5,1,2); plot(d{1}); title('D-1')
subplot(5,1,3); plot(d{2}); title('D-2')
subplot(5,1,4); plot(d{3}); title('D-3')
subplot(5,1,5); plot(d{4}); title('D-4')
subplot(6,1,6); plot(c_4); title('C-4')

%%%  FFFT Plot for Decomposition's Coeficents
Nfft=1024;
step=2*pi/Nfft;
t=-pi:step:pi-step;

figure;
for i=1:4
    temp=d{i};
    Ffft_d{i}=abs(fftshift(fft(temp,Nfft)));
    subplot(6,1,i);
    plot(t,Ffft_d{i});
    str = sprintf('d-%d ',i);
    title(str)
end
Ffft_c4=abs(fftshift(fft(c_4,Nfft)));
Ffft_f=abs(fftshift(fft(f,Nfft)));

subplot(615);
plot(t,Ffft_c4); title('c-4');
subplot(616);
plot(t,Ffft_f); title('f');



%%%% prob 4
d=[]; temp=[];
Ts=1/1000;
t=(0:Ts:2);

f=cos(2*pi*50*t);
p=(0.*(t<=0.5))+(0.25.*(t>0.5&t<=0.75))+(0.5.*(t<=1.25&t>0.75))+(0.25.*(t<=1.5&t>1.25))+(0.*(t>1.5&t<=2));
s=f+p;

[c,l]=wavedec(s,4,'db1');

d{4}=c(l(1)+1:l(1)+l(2));
d{3}=c(l(1)+l(2)+1:l(1)+l(2)+l(3));
d{2}=c(l(1)+l(2)+l(3)+1:l(1)+l(2)+l(3)+l(4));
d{1}=c(l(1)+l(2)+l(3)+l(4)+1:l(1)+l(2)+l(3)+l(4)+l(5));
c_4=c(1:l(1));


% %%%  FFFT Plot for Decomposition's Coeficents
Nfft=1024;
step=2*pi/Nfft;
t=-pi:step:pi-step;

figure;
for i=1:4
    temp=d{i};
    Ffft_d{i}=abs(fftshift(fft(temp,Nfft)));
    subplot(6,1,i);
    plot(t,Ffft_d{i});
    str = sprintf('d-%d ',i);
    title(str)
end
Ffft_c4=abs(fftshift(fft(c_4,Nfft)));
Ffft_f=abs(fftshift(fft(s,Nfft)));

subplot(615);
plot(t,Ffft_c4); title('c-4');
subplot(616);
plot(t,Ffft_f); title('s');


%%%% Plot Decomposition's Coeficents
figure;
subplot(6,1,1); plot(s); title('original signal')
subplot(6,1,2); plot(d{4}); title('D-4')
subplot(6,1,3); plot(d{3}); title('D-3')
subplot(6,1,4); plot(d{2}); title('D-2')
subplot(6,1,5); plot(d{1}); title('D-1')
subplot(6,1,6); plot(c_4); title('C_4')


%%%% recounstration and eliminate rectunguale pulse

c_4_New=zeros(1,length(c_4));

C=[c_4_New d{4} d{3} d{2} d{1}];
L=[length(c_4_New) length(d{4}) length(d{3}) length(d{2}) length(d{1}) length(f)];
f_New=waverec(C,L,'db1');



figure;
step=2/length(f);
tt=0:step:2-step;
subplot(211)
plot(tt,f_New);  title('reconstructed signal from p+f')
subplot(212);
plot(tt,s);      title(' signal  p+f')
