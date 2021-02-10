clc
clear all
close all
% f=zeros(100,100);
% f(1:100,1:49)=0;
% f(1:49,49:100)=0;
% f(50,50)=1;
% f(50,51:100)=2;
% f(51:100,50)=2;
% f(51:100,51:100)=3;

ff=imread('bloodcel_95.jpg');
ff=rgb2gray(ff);
ff=im2double(ff);
[M,N]=size(ff);
g=ff+0.1*randn(M,N);
figure
imshow(g);
w=4;
Gause=fspecial('gaussian',w,w/6);
f=imfilter(g,Gause);
figure,imshow(f,[]),title('filtered image');

[M,N]=size(f);
dx = [-1 0 1; -2 0 2; -1 0 1];
dy =dx';

dfx=imfilter(f,dx);
dfy=imfilter(f,dy);

%%%%%%%%%%%%%%%%% finding theta and gradian %%%%%%%%%%%%%
for m=1:M
    for n=1:N
        theta(m,n)=atan(dfy(m,n)/dfx(m,n));
        T(m,n)=(180/pi)*theta(m,n); 
        gradian(m,n)=sqrt((dfx(m,n))^2+(dfy(m,n))^2);
    end
end

%%%%%%%%%%%%%%%%%%%%%% theta quantization and finding thetahat %%%%%%%%%%%%%%%%%%%%%%%%%%%%
 for m=1:M
    for n=1:N
         if     T(m,n)<22.5,                  T1(m,n)=0;
         elseif T(m,n)>=22.5  & T(m,n)<67.5,  T1(m,n)=45;
         elseif T(m,n)>=67.5  & T(m,n)<112.5, T1(m,n)=90;
         elseif T(m,n)>=112.5 & T(m,n)<157.5, T1(m,n)=135;
         elseif T(m,n)>=157.5 & T(m,n)<180,   T1(m,n)=0;
         end
    end
 end
 
 %%%%%%%%%%%%%%%%%%%%%%% non-maximum supression %%%%%%%%%%%%%%%%%%%%%%%%%%%
 Mf=zeros(M,N);
 for m=2:M-1
    for n=2:N-1
            if T1(m,n)==0   & gradian(m,n)>=gradian(m-1,n)   & gradian(m,n)>=gradian(m+1,n),  Mf(m,n)=gradian(m,n);end
            if T1(m,n)==45  & gradian(m,n)>=gradian(m-1,n-1) & gradian(m,n)>=gradian(m+1,n+1),Mf(m,n)=gradian(m,n);end
            if T1(m,n)==90  & gradian(m,n)>=gradian(m,n-1)   & gradian(m,n)>=gradian(m,n+1),  Mf(m,n)=gradian(m,n);end 
            if T1(m,n)==135 & gradian(m,n)>=gradian(m-1,n+1) & gradian(m,n)>=gradian(m+1,n-1),Mf(m,n)=gradian(m,n);end
    end
 end

 %%%%%%%%%%%%%%%%%%%%%%%%% edge detection %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 R=zeros(M,N);
TH=0.8;
TL=0.7;
 for m=1:M
    for n=1:N
            if Mf(m,n)>=TH
                R(m,n)=1;
            end
    end
 end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 for m=2:M-1
    for n=2:N-1
            if T1(m,n)==0   & Mf(m-1,n)>TL   && Mf(m+1,n)>TL  ,  R(m,n)=1;end
            if T1(m,n)==45  & Mf(m-1,n-1)>TL & Mf(m+1,n+1)>TL , R(m,n)=1;end
            if T1(m,n)==90  & Mf(m,n-1)>TL   & Mf(m,n+1)>TL  ,  R(m,n)=1;end 
            if T1(m,n)==135 & Mf(m-1,n+1)>TL & Mf(m+1,n-1)>TL , R(m,n)=1;end
    end
 end
 figure,imshow(R),title(['canny ede detection with''TL=',num2str(TL),' TH=',num2str(TH)]);  
