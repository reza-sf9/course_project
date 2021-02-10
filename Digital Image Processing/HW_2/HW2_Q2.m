clear all
close all
clc
L=16;
b2=imread('baby2.jpg');
b2g=rgb2gray(b2);
figure
imshow(b2g,[])
title('baby2 original picture')

for i=0:255
   hist(i+1) = length(find(b2g == i));
end
figure,plot(hist,'LineWidth',2),title('Histogram of orginal Image')

%% uniform quantization
f_uniform=fix(b2g/(256/L))*16;
figure,imshow(f_uniform)
title(['uniform quantization with ',num2str(L),' level'])
r_uniform=0:16:240
for i=1:length(r_uniform)
    hist_uniform(i) = length(find( f_uniform == r_uniform(i) ));
end
figure, plot(r_uniform,hist_uniform,'b'),hold on
plot(r_uniform,hist_uniform,'k.','MarkerSize',30) %histogram of uniform quantization
title('histogram of picture after uniform Quamtization')
ti=fix([0:255/L:255]);
disp('ti of uniform quantization')
disp(ti)
[M,N]=size(f_uniform);
%% Max lloyd quantization
p=0;
h=imhist(b2g);
error=100;
while error>.1
    p=p+1;
    for i=1:L
        v=h(fix(ti(i):ti(i+1))+1);
        r(i)=[ti(i):ti(i+1)]*v/sum(v);
    end
    
    for j=2:L
        tt(j)=.5*(r(j-1)+r(j));
    end
    
    tt(L+1)=255;
    error=sum((tt-ti).^2);
    ti=tt;
    if error<.1
        break;
    end
    
end

f_LLoyd = zeros(M,N);
for l=2:L+1
    temp = find( ti(l-1)<b2g & b2g < ti(l));
    f_LLoyd(temp) = fix(r(l-1));
end
figure
imshow(f_LLoyd,[])
title(['LLoyd-max quantization with ',num2str(L),' level'])
fix_r = fix(r);
for i=1:length(fix_r)
    hist_LLoyd(i) = length(find( f_LLoyd == fix_r(i) ));
end
figure, plot(fix_r,hist_LLoyd,'r'),hold on
plot(fix_r,hist_LLoyd,'k.','MarkerSize',30),title('histogram of picture after maxlloyd quantization')
disp('ti of max lloyd quantization')
disp(fix(ti));
disp('ri of max lloyd quantization')
disp(r);
[M,N] = size(b2g);

f=zeros(M,N);
for i=0:255
   temp = find(b2g==i);
   f(temp) = i;
end 
f_uniform1 =zeros(M,N);
for i=0:16:240
   temp =  find(f_uniform==i);
   f_uniform1(temp) = i;
end



error_LLoyd =    ((norm(f(:)-f_LLoyd(:)))/norm(f(:)))*100;

error_Uniform =  ((norm(f(:)-f_uniform1(:)))/norm(f(:)))*100;