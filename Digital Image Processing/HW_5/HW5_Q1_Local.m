%%%%%%local Histogram Equalization%%%%%%%%%%
clc,clear all,close all;
f=imread('lowcontrast.jpg');
f=rgb2gray(f);
figure,imshow(f),title('original image');
[M,N]=size(f);
w=3;

% fhat=zeros(M+2*fix(w/2),N+2*fix(w/2));
% fhat(fix(w/2)+1:M1-fix(w/2),fix(w/2)+1:N1-fix(w/2))=f;

fhat = f;
[M1,N1]=size(fhat);
%**finding the histogram below the window
for m=fix(w/2)+1:M1-fix(w/2)-1
    for n=fix(w/2)+1:N1-fix(w/2)-1
        
        h=zeros(1,256);
        for i=-fix(w/2):fix(w/2)
            for j=-fix(w/2):fix(w/2)
                r=fhat(m+i,n+j);
                h(r+1)=h(r+1)+1;
            end
        end
        hist=h/sum(h);
        %****finding the cumulative function below the window
        F=zeros(1,256);
        for t=1:256
            F(t)=sum(hist(1:t));
        end
        s=fhat(m,n);
        g(m,n)=F(s+1);
    end
end

figure,imshow(g,[]),title(['equalized image by local histogram equalization,w= ',num2str(w)]);
