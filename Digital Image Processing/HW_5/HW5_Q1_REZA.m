clc;
clear;
close all;
%%%% Reza Eddition
%%%% Histogram Equalization (Global Method)
f=imread('lowcontrast.jpg');
f=rgb2gray(f);
f=im2double(f).*255;
figure,imshow(f,[]),title('Original Image')
[M,N] = size(f);
Pr_Dist=zeros(1,256);   %%% Probliility Distribution
for i=0:255
    Pr_Dist(1,i+1) = length(find(f==i))/(M*N);
end
% figure,plot(Pr_Dist),title('Probliility Distribution of Original Image');

Pr_Dens = zeros(1,256); %%% Probliility Density
for i=1:256
    Pr_Dens(i)=sum(Pr_Dist(1,1:i));
end
% figure,plot(Pr_Dens),title('Probliility Distribution of Original Image')

f_Global=zeros(M,N);
for i=0:255
    ind_temp = find(f==i);
    f_Global(ind_temp) = 255*Pr_Dens(i+1);
end
figure,imshow(f_Global,[]),title('Improved Contrast with Global Method');

figure,
subplot(211),histogram(f,256),title('Histogram of Original Image')
subplot(212),histogram(f_Global,256),title('Histogram after Equlaization')

%%%%% Local Histogram

for W=15:2:15 

w_2 = floor(W/2);

% f_Local = zeros(M-2*(w_2+1),N-2*(w_2+1));
p=0;
for m = w_2 +1 : M-w_2-1
    p=p+1; q=0;
    for n = w_2 +1 : N-w_2-1
        q=q+1;
        f_temp = f(m-w_2:m+w_2 , n-w_2:n+w_2);
        
        Pr_Dist_temp=zeros(1,256);   %%% Probliility Distribution
        for i=0:255
            Pr_Dist_temp(1,i+1) = length(find(f_temp==i))/(W^2);
        end
%         close all;
%         figure,stem(Pr_Dist_temp)
        
        Pr_Dens_temp = zeros(1,256); %%% Probliility Density
        for i=1:256
            Pr_Dens_temp(i)=sum(Pr_Dist_temp(1,1:i));
        end
%         figure,plot(Pr_Dens_temp)
        
       f_Local(p,q) = Pr_Dens_temp(1,1+f(m,n));
    end
end

gcf=figure;
imshow(f_Local),title(['Improved Image after Local Equalization - with Window Length = ',num2str(W)]);
label_name = sprintf('local_W_%d.jpg',W);
saveas(gcf,label_name)
end