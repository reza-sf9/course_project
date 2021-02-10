clc;
clear;
close all;

%%%%%%%%%%%%%%%%%%%%%%%%%% Problem 3
syms x;
%%%%% Calculate Scale

%%%% Mexican hat
Mex_hat_p2 = @(x) (1/sqrt(2*pi)*(1-x.^2).*exp(-x.^2/2)).^2; %% mexican hat power 2
% figure
% subplot(211),ezplot(Mex_hat_p2,[-4 4 -.1 .2]), title('Mexican hat power 2 , with \sigma = 1') , grid minor
int_psi_p2 = integral(Mex_hat_p2,-Inf,Inf);

%%%%% recived signal
s0=2; t0=1;
rec_ss = @(x) (1/sqrt(2*pi)*(1-(1/s0*(x-t0)).^2).*exp(-(1/s0*(x-t0)).^2/2)).^2; %% mexican hat shift 1; scale 2;  power 2

% subplot(212),ezplot(rec_ss,[-7  7 .2]), title('Recived Signal - Mexican hat power 2, with \sigma = 1 and scale factor = 3'), grid minor
int_rec_sig = integral(rec_ss,-Inf,Inf);

S = int_rec_sig/int_psi_p2;



%%%%% Calculate  Shift
Mex_hat_p2_s2 = @(x) (1/sqrt(2*pi)*(1-(x/2).^2).*exp(-(x/2).^2/2)).^2; %% mexican hat power 2 , scale 2 
int_psi_p2_s2 = integral(Mex_hat_p2_s2,-Inf,Inf);                      


int_rec_shiftPsi=0;
t0=-.01;
while abs(int_rec_shiftPsi - int_psi_p2_s2) > 0.00001  
 t0=t0+.01;
Mex_hat_s2 = @(x) 1/sqrt(2*pi)*(1-(1/S*(x-t0)).^2).*exp(-(1/S*(x-t0)).^2/2); %% mexican hat scale 2
% figure;
ezplot(Mex_hat_s2,[-10 10]);
%%%%% recived signal

rec_ss = @(x) 1/sqrt(2*pi)*(1-(1/2*(x-1)).^2).*exp(-(1/2*(x-1)).^2/2); %% mexican hat shift 1; scale 2;

hold on, ezplot(rec_ss,[-10 10]);
legend('mexican hat with scaling','recived signal');
title(sprintf(' shift factor = %f',t0));

%%% fun is Diffrencce f(t) * ?((t-?)/2)
fun = @(x) (1/sqrt(2*pi)*(1-(1/S*(x-t0)).^2).*exp(-(1/S*(x-t0)).^2/2)).*(1/sqrt(2*pi)*(1-(1/2*(x-1)).^2).*exp(-(1/2*(x-1)).^2/2));

int_rec_shiftPsi = integral(fun,-Inf,Inf);

end
k


