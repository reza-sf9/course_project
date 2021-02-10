clc
clear
close all

%%%%%%% Create and Image in movement 
r=imread('C:\Users\Reza_SF\Desktop\image\4.jpg');
r=rgb2gray(r);   
figure,imshow(r),title('Original Image')
f=im2double(r);  
[M,N] = size(f);

for a= 1:10
    for b=1:10
        T=40; 
        F=fft2(f);

        for v=1:N
            for u=1:M
                o= (a*u + b*v)/2 ;  
                H(u,v)= (T/o)   *  exp((-o) * 1i)  * sind(o);
                G1(u,v) = H(u,v) * F(u,v);
            end
        end
        
        %%%%% adding Noise
        n=0.003 * rand(M,N);
        Noise=fft2(n);
        
        G = G1 + Noise;  
        figure,imshow(real(ifft2(G))),title(['Noisy Image a= ',num2str(a),'  b= ',num2str(b)])
    end
end


% %%%%%% Problem Image
% load('restoration.mat');
% figure, imshow(g),title('Original Image')
% G = fft2(g);
% [M,N] = size(g);



%%%% Generate Laplacian  for CLS Method
Lap = zeros(M,N);
temp = ones(3,3);
temp(2,2) = -8;
Lap(1:3,1:3) =temp;
LAP = fft2(Lap);
LAP2 = abs(LAP).^2;


for a=2:2
    for b=13:1:13
        for T=50:10:50
            H=[];
            for u=1:M
                for v=1:N
                    coef = (a*u + b*v)/2;
                    H(u,v) = (T/coef) * exp(-1i*coef) * sind(coef);
                end
            end
            H_conj = conj(H);
            H2 = H.*H_conj;
            
            %%%%%%%% Wener Method
            
            for k_Wiener=0.01:.01:0.01
                Wiener_Filt = H_conj./(H2 + k_Wiener);
                
                F_hat_Winer = G.* Wiener_Filt;
                f_hat_Winer = real(ifft2(F_hat_Winer));
                
                str_tit = sprintf('Restored Image usding Wiener Filter \n Filter Parameters : a=%d b=%d k=%s T=%d',...
                    a,b,num2str(k_Wiener),T);
                figure,imshow(f_hat_Winer),title(str_tit)
            end
            
            %%%%%%% CLS Method
            for gamma=.001:.001:.001
                CLS_Filt = H_conj ./ (H2 + gamma.*(LAP2));
                F_hat_CLS = G.* CLS_Filt;
                f_hat_CLS = real(ifft2(F_hat_CLS));
                
                str_tit = sprintf('Restored Image usding CLS Filter \n Filter Parameters : a=%d b=%d \\gamma=%s',a,b,num2str(gamma));
                figure,imshow(f_hat_CLS),title(str_tit)
            end
            
        end
    end
end


