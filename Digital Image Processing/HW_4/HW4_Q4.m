clc;
clear;
close all;

%%%%%%%%%%%%%%%%%
f= imread('baby2.jpg');
f= rgb2gray(f);
f = im2double(f);
[M,N] = size(f);

figure, imshow(f,[]),title('Original Images');

%%% zeropad
a_r = mod(M,8); a_c = mod(N,8);
if a_r ~= 0
   a_r = 8 - a_r; 
end
if a_c ~= 0
    a_c = 8 - a_c; 
end
temp = zeros(M+  a_r , N+a_c);
temp(1:M,1:N) = f;
ff= temp;
[M1,N1] = size(ff);


%%%%%%
f_4 = zeros(M1,N1);    f_2 = zeros(M1,N1);

for m=1:8:M1
    for n=1:8:N1
        DCT_4 = zeros(8,8); DCT_2 = zeros(8,8);
        
        DCT_temp = dct2(ff(m:m+7,n:n+7));
        DCT_4(1:4,1:4) =  DCT_temp(1:4,1:4); %%%% preserve only 4*4 left up 
        DCT_2(1:2,1:2) =  DCT_temp(1:2,1:2); %%%% preserve only 2*2 left up 
        
        f_4(m:m+7,n:n+7) = idct2(DCT_4);     %%%% reconstruct f
        f_2(m:m+7,n:n+7) = idct2(DCT_2);
    end
end

%%%%%% calulation normlized matrix  with padding
norm_ff = ff./norm(ff);
norm_f_2_pad = f_2./norm(f_2);
norm_f_4_pad  = f_4./norm(f_4);
%%%%%% caculation error rate
e_2_pad = norm(norm_ff - norm_f_2_pad )*100;
e_4_pad = norm(norm_ff - norm_f_4_pad )*100;

%%%%%% getting rid of padding
f_4 = f_4(1:M,1:N);
f_2 = f_2(1:M,1:N);
%%%%%% calulation normlized matrix  without padding
norm_f = f./norm(f);
norm_f_2 = f_2./norm(f_2);
norm_f_4 = f_4./norm(f_4);
%%%%%% caculation error rate
e_2 = norm(norm_f - norm_f_2)*100;
e_4 = norm(norm_f - norm_f_4)*100;


%%%%%% calulation compression rate
compression_rate_2 = 100*(1-(2^2/8^2));
compression_rate_4 = 100*(1-(4^2/8^2));
%%%%%% plotting Images
figure, imshow(f_4,[]), title(['using only 4*4 left up of DCT coefficients , compression rate is ',num2str(compression_rate_4),' %']);
figure, imshow(f_2,[]), title(['using only 2*2 left up of DCT coefficients, compression rate is ',num2str(compression_rate_2),' %']);