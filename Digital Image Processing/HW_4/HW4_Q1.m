clc;
clear;
close all;


load('exersize4.mat');

figure;
imshow(F)

[M,N] = size(F);
p=0;
for m=1:8:M
    p=p+1; q=0;
    for n=1:8:N
        q=q+1;
   temp = F(m:m+7,n:n+7); 
   DCT_temp = dct2(temp);
   index_41(p,q) = DCT_temp(4,1);        
   index_23(p,q) = DCT_temp(2,3);   
   index_32(p,q) = DCT_temp(3,2);   
    end
end

Message_1 = zeros(M/8,N/8);
Message_2 = zeros(M/8,N/8);

temp_1 = index_41 - index_23;
a_1=find(temp_1>=0);
Message_1(a_1) = 1;

temp_2 = index_41 - index_32;
a_2=find(temp_2>=0);
Message_2(a_2) = 1;

figure;
subplot(211),imshow(Message_1),title('if we choses pixel(2,3)');
subplot(212),imshow(Message_2),title('if we choses pixel(3,2)');

