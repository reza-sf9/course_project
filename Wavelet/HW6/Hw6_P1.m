clc;
clear;
close all;

%%%%%%%%%%%%%%%%%%%%%% Problem 1  %%%%%%%%%%%%%%%%%%%%%%%

%%% calculate coeficients of filters

[h g] = wfilters('db2','r');
h
figure;
subplot(211); stem(h);  title('coefficents of h filter (db4)')
subplot(212); stem(g);  title('coefficents of g filter (db4)')

%%% plot scaling function and wavelet

[phi,psi,xval] = wavefun('db4',10); %% for 10 iteration
figure;
subplot(211);
plot(xval,phi);
title('db4 Scaling Function (\phi)');
subplot(212);
plot(xval,psi);
title('db4 Wavelet (\psi)');