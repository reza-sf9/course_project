clc;
clear;
close all;



%%%%%% Prob1 %%%%%
wname='db2';

a=imread('cat.jpg');
a=rgb2gray(a);
a=im2double(a);

a = (a-min(a(:))) ./ (max(a(:)-min(a(:)))); %%%%    SCALED A BETWEEN 0-1

%%% Decomposition
[C,S] = wavedec2(a,1,wname);

cA=appcoef2(C,S,wname,1);            %% get Approximation Coeffitents(LL)

[cH,cV,cD] = detcoef2('all',C,S,1);  %% get Horizental Vertical & Diogonal Coefficents


% cA = (cA-min(cA(:))) ./ (max(cA(:)-min(cA(:)))); %%%%    SCALED A BETWEEN 0-1
% 
% 
% %%%% PLOTTING COMPONENTS
% figure;
% subplot(321)
% imshow(cA); title('LL');
% subplot(322)
% imshow(cH); title('LH');
% subplot(323)
% imshow(cV); title('HL');
% subplot(324)
% imshow(cD); title('HH');
% subplot(3,2,[5,6])
% imshow(a); title('pic');






%%%% calculate mean and std and determin threshold
k=.8;
mean_cH=mean(cH(:)); std_cH=std(cH(:)); 
Th_cH= mean_cH + k * std_cH;
cH_S=wthresh(cH,'s',Th_cH);     %% soft Thresholding
cH_H=wthresh(cH,'h',Th_cH);     %% hard Thresholding


mean_cV=mean(cV(:)); std_cV=std(cV(:)); 
Th_cV= mean_cV + k * std_cV;
cV_S=wthresh(cV,'s',Th_cV);     %% soft Thresholding
cV_H=wthresh(cV,'h',Th_cV);     %% hard Thresholding

mean_cD=mean(cD(:)); std_cD=std(cD(:)); 
Th_cD= mean_cD + k * std_cD;
cD_S=wthresh(cD,'s',Th_cD);     %% soft Thresholding
cD_H=wthresh(cD,'h',Th_cD);     %% hard Thresholding

Density_Ratio = (length(find(cH_S(:)==0)) + length(find(cV_S(:)==0)) + length(find(cD_S(:)==0)))/(3*length(cH(:)));


ST_a = idwt2(cA,cH_S,cV_S,cD_S,wname);   %%% reconstructed with Soft Threshold's coesfficents
HT_a = idwt2(cA,cH_H,cV_H,cD_H,wname);   %%% reconstructed with Hard Threshold's coesfficents

%%% SHOW IMAGES
figure;
subplot(131), imshow(a),     title('original image')
subplot(132), imshow(ST_a),  title('soft threshold')
subplot(133), imshow(HT_a),  title('Hard threshold')

%%% Normalize
norm_a=norm(a);
norm_ST_a=norm(ST_a);
norm_HT_a=norm(HT_a);

a_n=a/sqrt(norm_a);
ST_a_n=ST_a/sqrt(norm_ST_a);
HT_a_n=HT_a/sqrt(norm_HT_a);

e_ST=norm(a_n-ST_a_n);     %%% error of soft Thresholding
e_HT=norm(a_n-HT_a_n);     %%% error of hard Thresholding


%%%% PLOTTING COMPONENTS
% figure;
% subplot(321)
% imshow(cA); title('LL');
% subplot(322)
% imshow(cH); title('LH');
% subplot(323)
% imshow(cV); title('HL');
% subplot(324)
% imshow(cD); title('HH');
% subplot(3,2,[5,6])
% imshow(a); title('pic');


