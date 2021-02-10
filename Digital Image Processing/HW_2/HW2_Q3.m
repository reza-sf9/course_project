clc;
clear;
close all;

I = imread('text_10.jpg');
I = rgb2gray(I);
figure, imshow(I), title('Original Image');

hist = zeros(1,256);
for i=0:255
    hist(i+1) = length(find(I(:) == i));
end
figure, plot(0:255,hist, 'LineWidth',2);
title('Histogram of Original Image')

%%%%%%%%%%%%%%%% Binary, using T method
error = 1;
p = 0;
T=127;
while error > .01
    p= p +1;
    
    m_1(p) = (0:fix(T)-1)*(hist(1:fix(T))/sum(hist(1:fix(T)))).';
    m_2(p) = (fix(T):255)*(hist(fix(T)+1:256)/sum(hist(fix(T)+1:256))).';
    TT = .5*(m_1(p)+m_2(p));
    error(p) = norm(T-TT);
    T =TT;
end
T_T =T;
a= find(I>=0 & I<=T_T);  %%% indices that lower than T
b= find(I>T_T);          %%% indices taht greater than T
I_Binary_T  = I;
I_Binary_T(a) =0;
I_Binary_T(b) =255;
figure, imshow(I_Binary_T);
title('Binary Image using T  mehtod ')


%%%%%%%%%%% Mohammad Motamedi OTSU
%% OTSU
%  f= I; h = hist;
%  [M,N] =size(f);
% g1=zeros(M,N);
% [M,N]=size(f);
% L=length(h);
% for T=1:length(h)
%     sigma1=0;
%     sigma2=0;
%     T1=h(1:T);
%     T2=h(T+1:L);
%     sum1=sum(T1);
%     sum2=sum(T2);
%     p1=T1/sum1;
%     p2=T2/sum2;
%     mean1=sum(p1.*T1);
%     mean2=sum(p2.*T2);
%     for i=1:T
%         sigma1=sigma1+(p1(i)*((i-mean1).^2));    
%     end
%     for i=T+1:length(p2) %%%% why untill L?????
%         sigma2=sigma2+(p2(i)*((i-mean2).^2));
%     end
%     P1=sum1/sum(h);
%     P2=sum2/sum(h);
%     J(1,T)=P1*sigma1+P2*sigma2;
% end 
% figure(4);
% plot(J)
% newT=min(find(J==min(J)));
% for m=1:M
%     for n=1:N
%         if f(m,n)>=newT
%             g1(m,n)=1;
%         end
%     end
% end
% figure(5);
% imshow(g1,[]);title('OTSU');



%%%%%%%%%%%%%%%% Binary, using OTSU method
p = hist/sum(hist);
sum(hist)
for i=1:256
    T=i;
    q_1 = sum(p(1:T));
    q_2 = sum(p(T+1:256));
    
    p_1 = p(1:T)/q_1; p_2 = p(T+1:256)/q_2; 
    
    m_1 = (1:T)*p_1.';
%     m_1 = m_1;
    
    m_2 = (T+1:256)*p_2.';
%     m_2 = m_2;
    
    
    sigma_1 = ((((1:T)-m_1).^2)*p_1.');
    
    sigma_2 = (((T+1:256)-m_2).^2 * p_2.');
    
%     P_1 = q_1/sum(p);
%     P_2 = q_2/sum(p);
    
    J(i)=q_1* sigma_1 + q_2*sigma_2; 

end
min_index = find(J == min(J));
T_J = min_index(1);


figure, plot(J), xlim([0 256]),hold on;
plot(T_J,J(T_J),'r.','MarkerSize',30);
title('J Criteria')

a= find(I>=0 & I<=T_J);  %%% indices that lower than T
b= find(I>T_J);          %%% indices taht greater than T
I_Binary_OTSU  = I;
I_Binary_OTSU(a) =0;
I_Binary_OTSU(b) =255;

figure, imshow(I_Binary_OTSU);
title('Binary Image Using OTSU Method ')

level = graythresh(I);
mat_otsu = level*256;

%%%%%%%%%%%%%%%%%%%%%%%%%    Devide Image to 4 part  

[M,N] = size(I);
I_4{1} = I(1:fix(M/2),1:fix(N/2));           %%%% Top left Image
I_4{2} = I(1:fix(M/2),fix(N/2)+1:N);         %%%% Top right Image
I_4{3} = I(fix(M/2)+1:end,1:fix(N/2));       %%%% Bottom left Image
I_4{4} = I( fix(M/2)+1:end,fix(N/2)+1:N);    %%%% Bottom right Image


for i=1:4
    I_temp = I_4{i};
    for j=0:255
        temp(j+1) = length(find(I_temp(:) == j));
    end
    hist_4{i} = temp;
end
%%%%%%%%%%%%%%%%%%%%%%

%%%%% plotting Image after deviding to 4 part
figure,subplot(221),imshow(I_4{1}),title('Image - Top Left')
subplot(222),imshow(I_4{2}),title('Image - Top Right')
subplot(223),imshow(I_4{3}),title('Image - Bottom Left')
subplot(224),imshow(I_4{4}),title('Image - Bottom Right')

%%%%% Histogram of Image after deviding to 4 part

figure,subplot(221),plot(hist_4{1},'LineWidth',2),title('Histogram - Top Left'),xlim([1 256])
subplot(222),plot(hist_4{2},'LineWidth',2),title('Histogram - Top Right'),xlim([1 256])
subplot(223),plot(hist_4{3},'LineWidth',2),title('Histogram - Bottom Left'),xlim([1 256])
subplot(224),plot(hist_4{4},'LineWidth',2),title('Histogram - Bottom Right'),xlim([1 256])

%%%%%%%%%%%%%%%%% Binary metod 4 Image
error = 1;

for i=1:4
    hist=[];
    hist = hist_4{i}; 
    p = 1;
    T=127;
    e=[]; e=10;
    while e > .01        
      
        m_1 = (0:fix(T)-1)*(hist(1:fix(T))/sum(hist(1:fix(T)))).';
        m_2 = (fix(T):255)*(hist(fix(T)+1:256)/sum(hist(fix(T)+1:256))).';
        TT = .5*(m_1 + m_2);
        e(p) = norm(T-TT);
        T =TT;
        p= p +1;
    end
    T_T_4(i) = T;
    error_4{i} = e;
end

for i=1:4
    temp = I_4{i};
    a= find(temp>=0 & temp<=T_T_4(i));  %%% indices that lower than T
    b= find(temp > T_T_4(i));          %%% indices taht greater than T
    temp(a) =0;
    temp(b) =255;
    I_Binary_T_4{i}  = temp;
end


figure, subplot(221),imshow(I_Binary_T_4{1});title('Image Binary - Top Left')
subplot(222),imshow(I_Binary_T_4{2});title('Image Binary - Top Right')
subplot(223),imshow(I_Binary_T_4{3});title('Image Binary - Bottom Left')
subplot(224),imshow(I_Binary_T_4{4});title('Image Binary - Bottom Right')

I_B_4_T = [I_Binary_T_4{1} I_Binary_T_4{2}; I_Binary_T_4{3} I_Binary_T_4{4}]; %%% merge 4 Image together
figure , imshow(I_B_4_T),title('Merge 4 Image together anfter using Binary Method(T)')

%%%%%%%%%%%%%%%% OTSU  4 Image

for j=1:4
    hist = hist_4{j};
    p=hist/sum(hist);
    for i=1:256
    T=i;
    q_1 = sum(p(1:T));
    q_2 = sum(p(T+1:256));
    
    m_1 = (0:T-1)*p(1:T).';
    m_1 = m_1/q_1;
    
    m_2 = (T:255)*p(T+1:256).';
    m_2 = m_2/q_2;
    
    sigma_1 =0;
    sigma_1 = ((((0:T-1)-m_1).^2)*p(1:T).')/q_1;
    
    sigma_2 = 0;
    sigma_2 = (((T:255)-m_2).^2 * p(T+1:256).')/q_2;
    
    J_4(j,i)=q_1* sigma_1 + q_2*sigma_2; 

    end
end

for i=1:4
    min_index = find(J_4(i,:) == min(J_4(i,:)));
    T_J_4(i)= min_index(1);
end


figure,subplot(221),plot(J_4(1,:), 'LineWidth',2),title('J Criteria - Top Left'),hold on;
plot(T_J_4(1),J_4(1,T_J_4(1)),'r.','MarkerSize',30);
subplot(222),plot(J_4(2,:), 'LineWidth',2),title('J Criteria - Top Right'),hold on;
plot(T_J_4(2),J_4(2,T_J_4(2)),'r.','MarkerSize',30);
subplot(223),plot(J_4(3,:), 'LineWidth',2),title('J Criteria - Bottom Left'), hold on;
plot(T_J_4(3),J_4(3,T_J_4(3)),'r.','MarkerSize',30);
subplot(224),plot(J_4(4,:), 'LineWidth',2),title('J Criteria - Bottom Right'),hold on
plot(T_J_4(4),J_4(4,T_J_4(4)),'r.','MarkerSize',30);


for i=1:4
    a=[];b=[];
    a= find(I_4{i} >= 0 & I_4{i} <= T_J_4(i));   %%% indices that lower than T
    b= find(I_4{i} > T_J_4(i));                  %%% indices taht greater than T
    temp = I_4{i};
    temp(a) =0;
    temp(b) =255;
    I_Binary_J_4{i} = temp;
end

figure,subplot(221),imshow(I_Binary_J_4{1}),title('Image Binary(OTSU method) - Top Left')
subplot(222),imshow(I_Binary_J_4{2}),title('Image Binary(OTSU method) - Top Right')
subplot(223),imshow(I_Binary_J_4{3}),title('Image Binary(OTSU method) - Bottom Left')
subplot(224),imshow(I_Binary_J_4{4}),title('Image Binary(OTSU method) - Bottom Right')

I_B_J_4 = [I_Binary_J_4{1} I_Binary_J_4{2}; I_Binary_J_4{3} I_Binary_J_4{4}];
figure,imshow(I_B_J_4),title('Merge 4 Image together after using Binary Method(OTSU)')
