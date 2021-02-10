close all
clear
clc

g=imread('finger.tif');
%%%% convert image 2 square image
[M , N] = size(g);
temp = [M N];
ind_min = find(temp == min(temp));
switch ind_min
    case 1
        temp1 = N-M;
        g1 = g(1:end , 1+temp1/2 : end - temp1/2);
    case 2
        temp1 = M-N;
        g1 = g(1+temp1/2 : end - temp1/2 , 1:end);
end

g1=im2double(g1);
figure
imshow(g1,[])
title('Orginal Image (square size)')

%% Generate QUINCUNX 8 filter
w=size(g1,1);
s=w/2;
h1=zeros(s);

figure
square1=zeros(w);
square1(1:s-1,1:s-1)=ones(s-1);
square1(s+2:w,s+2:w)=ones(s-1);subplot(221),imshow(square1,[]),title('square1')
square2=zeros(w);
square2(1:s-1,s+2:w)=ones(s-1);
square2(s+2:w,1:s-1)=ones(s-1);subplot(222),imshow(square2,[]),title('square2')
h2=h1;
for m=1:s
    for n=1:s
        if m<n
            h1(m,n)=1;
            
        end
        if m>s-n
            h2(m,n)=1;
        end
    end
end

f1=[h1,1-h2;h2,1-h1];  subplot(223),imshow(f1,[]),title('f1')
f2=[1-h1,h2;1-h2,h1];  subplot(224),imshow(f2,[]),title('f2')


h1=zeros(w);
h2=h1;
for m=1:w
    for n=1:w
        if m<n
            h1(m,n)=1;
            
        end
        if m>w-n
            h2(m,n)=1;
        end
    end
end

rect1=[h1,1-h2;h2,1-h1];
rect2=[1-h1,h2;1-h2,h1];

for m=-s:s-1
    for n=-s:s-1
        y11(m+s+1,n+s+1)=rect1(m+n+w+1,n+w+1);
        y12(m+s+1,n+s+1)=rect1(m-n+w+1,n+w+1);
        y13(m+s+1,n+s+1)=rect1(m+w+1,n+m+w+1);
        y14(m+s+1,n+s+1)=rect1(m+w+1,n-m+w+1);
    end
end
figure,subplot(221),imshow(y11,[]),title('y11')
subplot(222),imshow(y12,[]),title('y12')
subplot(223),imshow(y13,[]),title('y13')
subplot(224),imshow(y14,[]),title('y14')

for m=-s:s-1
    for n=-s:s-1
        y21(m+s+1,n+s+1)=rect2(m+n+w+1,n+  w+1);
        y22(m+s+1,n+s+1)=rect2(m-n+w+1,n+  w+1);
        y23(m+s+1,n+s+1)=rect2(m+  w+1,n+m+w+1);
        y24(m+s+1,n+s+1)=rect2(m+  w+1,n-m+w+1);
    end
end

figure,subplot(221),imshow(y21,[]),title('y21')
subplot(222),imshow(y22,[]),title('y22')
subplot(223),imshow(y23,[]),title('y23')
subplot(224),imshow(y24,[]),title('y24')

figure,
z1=f1.*square1;  subplot(221),imshow(z1,[]),title('z1')
z2=f1.*square2;  subplot(222),imshow(z2,[]),title('z2')
z3=f2.*square1;  subplot(223),imshow(z3,[]),title('z3')
z4=f2.*square2;  subplot(224),imshow(z4,[]),title('z4')


figure
Q_8(:,:,1)=y11.*square2;   subplot(421),imshow(Q_8(:,:,1),[]),title('q1')
Q_8(:,:,2)=y12.*square1;   subplot(422),imshow(Q_8(:,:,2),[]),title('q2')
Q_8(:,:,3)=y13.*z4;        subplot(423),imshow(Q_8(:,:,3),[]),title('q3')
Q_8(:,:,4)=y14.*z3;        subplot(424),imshow(Q_8(:,:,4),[]),title('q4')
Q_8(:,:,5)=y21.*z2;        subplot(425),imshow(Q_8(:,:,5),[]),title('q5')
Q_8(:,:,6)=y22.*z1;        subplot(426),imshow(Q_8(:,:,6),[]),title('q6')
Q_8(:,:,7)=y23.*square2;   subplot(427),imshow(Q_8(:,:,7),[]),title('q7')
Q_8(:,:,8)=y24.*square1;   subplot(428),imshow(Q_8(:,:,8),[]),title('q8')
suptitle('8 quincunx filter');


%%% this part causes the frequncy response be symmetric around 0
for m=0:w-1
    for n=0:w-1
        f(m+1,n+1)=(-1)^(m+n)*g1(m+1,n+1);
    end
end
F=fft2(f);
for k=1:8
    Fhat(:,:,k)=Q_8(:,:,k).*F;
end
[M , N] = size(Q_8(:,:,1));

%%%% plot each part of 16 part of Image after filtering
for k=1:8
    fhat=ifft2(Fhat(:,:,k));
    fhat=real(fhat);
    for m=0:w-1
        for n=0:w-1
            final_f(m+1,n+1)=(-1)^(m+n)*fhat(m+1,n+1);
        end
    end
    final_ff(:,:,k) = final_f;
%     figure,imshow(final_f,[]),title([num2str(k),'th part of Image after using QUINQUNX (8 Filter)'])
end

% 
% %%%%%% Fusing Image with using BLOCK Method 
% W = 5;
% w_2 = floor(W/2);
% 
% for k=1 : 8
%     std_f = zeros(1, 8);
%     for m= w_2 +1 : W : M - w_2 -1
%         for n= w_2 +1 : W : N - w_2 -1
%             for k=1 : 8
%                 temp = final_ff(m-w_2 : m+w_2 , n-w_2 : n+w_2,k);
%                 std_f(1,k) = std2(temp);
%             end
%             ind_max = min(find( std_f == max(std_f)));
%             f_fuse_block_8(m-w_2 : m+w_2 , n-w_2 : n+w_2) = final_ff(m-w_2 : m+w_2 , n-w_2 : n+w_2 ,ind_max);
%         end
%     end
% end
% 
% figure,imshow(f_fuse_block_8,[]),title('Fusing Image Using Block Method-(QUINCUNX 8 FILTER)')
% 
% 
% %%%%%% Fusing Image with using PIXEL Method 
% W = 5;
% w_2 = floor(W/2);
% 
% for k=1 : 8
%     std_f = zeros(1, 8);
%     for m= w_2 +1 : M - w_2 -1
%         for n= w_2 +1 : N - w_2 -1
%             for k=1 : 8
%                 temp = final_ff(m-w_2 : m+w_2 , n-w_2 : n+w_2,k);
%                 std_f(1,k) = std2(temp);
%             end
%             ind_max = min(find( std_f == max(std_f)));
%             f_fuse_pix_8(m , n) = final_ff(m , n ,ind_max);
%         end
%     end
% end
% 
% figure,imshow(f_fuse_pix_8,[]),title('Fusing Image Using Pixel Method-(QUINCUNX 8 FILTER)')
% 


%% %%%%%%%%%%%%%%%%%%%%%%%%% generate QUINCUNX 8 filter
[M,N] = size(Q_8(:,:,1));
p=0;
for i=1:2:16
    p=p+1;
    temp_m = Q_8(1,:,p);
    ind_first_1_m = min(find(temp_m==1));
    temp_n = Q_8(:,1,p);
    ind_first_1_n = min(find(temp_n==1));
    
    if  isempty(ind_first_1_m)==1
        ind_first_1_m = 480;
    elseif isempty(ind_first_1_n)==1
        ind_first_1_n = 480;
    end
    
    if ind_first_1_m < ind_first_1_n
        type_m_n = 'm';
    else
        type_m_n = 'n';
    end
    
    switch type_m_n
        
        case 'm'
            
            for m=1:M/2
                temp = Q_8(m,:,p);
                ind_first_1 = min(find(temp==1));
                ind_last_1 = max(find(temp==1));
                ind_med = floor(.5*(ind_first_1+ind_last_1));
                Q_16(m,ind_first_1:ind_med,i) = 1;
                Q_16(m,ind_med+1:ind_last_1,i+1) = 1;
            end
            
            for m=M/2+1:M
                temp = Q_8(m,:,p);
                ind_first_1 = min(find(temp==1));
                ind_last_1 = max(find(temp==1));
                ind_med = floor(.5*(ind_first_1+ind_last_1));
                Q_16(m,ind_first_1:ind_med,i+1) = 1;
                Q_16(m,ind_med+1:ind_last_1,i) = 1;
            end
            
        case 'n'
            for n=1:N/2
                temp = Q_8(:,n,p);
                ind_first_1 = min(find(temp==1));
                ind_last_1 = max(find(temp==1));
                ind_med = floor(.5*(ind_first_1+ind_last_1));
                Q_16(ind_first_1:ind_med,n,i) = 1;
                Q_16(ind_med+1:ind_last_1,n,i+1) = 1;
            end
            
            for n=N/2+1:N
                temp = Q_8(:,n,p);
                ind_first_1 = min(find(temp==1));
                ind_last_1 = max(find(temp==1));
                ind_med = floor(.5*(ind_first_1+ind_last_1));
                Q_16(ind_first_1:ind_med,n,i+1) = 1;
                Q_16(ind_med+1:ind_last_1,n,i) = 1;
            end
    end
end

%%%%%% Plot 16 QUINCUNX Filter
figure;
for i=1:2:8
    subplot(4,2,i),imshow(Q_16(:,:,i)),title(['q',num2str(i)]);
    subplot(4,2,i+1),imshow(Q_16(:,:,i+1)),title(['q',num2str(i+1)]);
end
suptitle('1-8 filters of 16 filetr of QUINCUNX Filter Bank');
figure;
for i=9:2:16
    subplot(4,2,i-8),imshow(Q_16(:,:,i)),title(['q',num2str(i)]);
    subplot(4,2,i-8+1),imshow(Q_16(:,:,i+1)),title(['q',num2str(i+1)]);
end
suptitle('9-16 filters of 16 filetr of QUINCUNX Filter Bank');


%%%%% imfilter
for k=1:16
    Fhat(:,:,k)=Q_16(:,:,k).*F;
end
[M , N] = size(Q_16(:,:,1));

%%%% plot each part of 16 part of Image after filtering
for k=1:16
    fhat=ifft2(Fhat(:,:,k));
    fhat=real(fhat);
    for m=0:w-1
        for n=0:w-1
            final_f(m+1,n+1)=(-1)^(m+n)*fhat(m+1,n+1);
        end
    end
    final_ff(:,:,k) = final_f;
    gcf,figure,imshow(final_f,[]);figure,imshow(final_f,[]),title([num2str(k),'th part of Image after using QUINQUNX (16 Filter)'])
end


%%%%%% Fusing Image with using BLOCK Method 
W = 5;
w_2 = floor(W/2);

for k=1 : 16
    std_f = zeros(1, 16);
    for m= w_2 +1 : W : M - w_2 -1
        for n= w_2 +1 : W : N - w_2 -1
            for k=1 : 16
                temp = final_ff(m-w_2 : m+w_2 , n-w_2 : n+w_2,k);
                std_f(1,k) = std2(temp);
            end
            ind_max = min(find( std_f == max(std_f)));
            f_fuse_block_16(m-w_2 : m+w_2 , n-w_2 : n+w_2) = final_ff(m-w_2 : m+w_2 , n-w_2 : n+w_2 ,ind_max);
        end
    end
end

figure,imshow(f_fuse_block_16,[]),title('Fusing Image Using Block Method-(QUINCUNX 16 FILTER)')


%%%%%% Fusing Image with using PIXEL Method 
W = 5;
w_2 = floor(W/2);

for k=1 : 16
    std_f = zeros(1, 16);
    for m= w_2 +1 : M - w_2 -1
        for n= w_2 +1 : N - w_2 -1
            for k=1 : 16
                temp = final_ff(m-w_2 : m+w_2 , n-w_2 : n+w_2,k);
                std_f(1,k) = std2(temp);
            end
            ind_max = min(find( std_f == max(std_f)));
            f_fuse_pix_16(m , n) = final_ff(m , n ,ind_max);
        end
    end
end

figure,imshow(f_fuse_pix_16,[]),title('Fusing Image Using Pixel Method-(QUINCUNX 16 FILTER)')
