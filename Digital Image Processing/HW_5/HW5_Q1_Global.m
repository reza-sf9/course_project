clc
clear all
close all
f=imread('lowcontrast.jpg');
f=rgb2gray(f);
figure,imshow(f),title('original image');
[M,N]=size(f);
%% %%%%%%%%%%%structing the histogam function%%%%%%%%%%%
h=zeros(1,256);
for m=1:M
    for n=1:N
        r=f(m,n);
        h(r+1)=h(r+1)+1;%h is the histogram function
    end
end
hist=h/sum(h);
figure
plot(hist)
title('histogram');
%%%%%%%%%structing the cumulative function%%%%%%%%%5
F=zeros(1,256);
for k=1:256
    F(k)=sum(hist(1:k));
end
figure,plot(F),title('cumulative function');
%%%%%%%%%%%implementing Histogram equalizeion %%%%%%%%%%55
g=zeros(M,N);
for m=1:M
    for n=1:N
        s=f(m,n);
        g(m,n)=F(s);%g is the resulted image after histogram equalization
    end
end
figure,imshow(g),title('equalized image by Global Histogram Equalization');

