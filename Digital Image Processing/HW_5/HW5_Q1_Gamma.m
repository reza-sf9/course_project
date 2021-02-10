clc
clear all
close all
f=imread('lowcontrast3.jpg');
f=rgb2gray(f);
f=im2double(f);
figure,imshow(f),title('original image');
[M,N]=size(f);
h=zeros(1,256);
for m=1:M
    for n=1:N
        r=f(m,n)*255;
        h(r+1)=h(r+1)+1;
    end
end
hist=h/sum(h);
%figure,plot(hist),title('histogram');
%%%%%%%%%%%%% gma correction %%%%%%%%%%
gama=4;
R=zeros(1,256);
S=zeros(1,256);
g=zeros(M,N);
for m=1:M
    for n=1:N
        r=f(m,n)*255;
        R(r+1)=r;
        s=(255^(1-gama))*(r^gama);
        S(r+1)=s;
        g(m,n)=s;
    end
end

figure,imshow(g,[]),title(['equalized image by gama correction,gama= ',num2str(gama)]);
