clc
clear;
close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%     problem 4

wname='db10';

step=2/500;
t=0:step:2-step;

%%%% generate signal
f = cos(2*pi*10*t);       

%%%%  generate rectangular noise
t0=1/20;
s=1;
p=zeros(1,500);
for i=0:step:2-step
    if mod(floor(s/20),2)==1
        p(s)=.2;
    else
        p(s)=-.2;
    end
    s=s+1;
end


Lambda_d1=zeros(1,10); Lambda_d2=zeros(1,10); Lambda_d3=zeros(1,10);
Lambda_d4=zeros(1,10); Lambda_d5=zeros(1,10); Lambda_d6=zeros(1,10);
Lambda_c_d=zeros(1,10);

%%% we use this FOR , for consider variation of Gussian Noise
for i=1:10
    
    mu = 0; std=.4;
    v=random('norm',mu,std,1,500); %% generate Gussian noise
    
    x = f + p + v;    %%%  signal + noise 
    
    %%%%%%%%%%% Decomposition signal
    lev=4;
    [c,l]=wavedec(x, lev ,wname);
    
    [d1 ,d2, d3, d4]=detcoef(c,l,[1 2 3 4 ]);
    [c_d]=appcoef(c,l,wname,lev);
    
    
    
    %%%% calculate Lambda for each level
    Lambda_d1(i) = calc_Lambda_Optimum(d1,'d1');
    Lambda_d2(i) = calc_Lambda_Optimum(d2,'d2');
    Lambda_d3(i) = calc_Lambda_Optimum(d3,'d3');
    Lambda_d4(i) = calc_Lambda_Optimum(d4,'d4');
    Lambda_c_d(i) = calc_Lambda_Optimum(c_d,'c_d');
end

th_type='s';
Lambda_d1=mean(Lambda_d1);                %% average of Lambda
d1_N = wthresh(d1,th_type,Lambda_d1);     %%  Thresholding

Lambda_d2=mean(Lambda_d2);
d2_N = wthresh(d2,th_type,Lambda_d2);     %%  Thresholding

Lambda_d3=mean(Lambda_d3);
d3_N = wthresh(d3,th_type,Lambda_d3);     %%  Thresholding

Lambda_d4=mean(Lambda_d4);
d4_N = wthresh(d4,th_type,Lambda_d4);     %%  Thresholding


Lambda_c_d =mean(Lambda_c_d);
c_d_N = wthresh(c_d,th_type,Lambda_c_d);      %%  Thresholding
 

%%%%% Reconstruction Signal after thresholding components 
% % %%% 4 level
C=[c_d_N d4_N d3_N d2_N d1_N];
L=[length(c_d_N) length(d4_N) length(d3_N) length(d2_N) length(d1_N) length(x)];
x_New=waverec(C,L,wname);



%%% plottong result of recondtruction
figure;
plot(t,x), title('signal + noise');
figure;
plot(t,x_New), title(sprintf('reconstructed signal , wavelet type = %s * Number of Level = %d',wname,lev));

