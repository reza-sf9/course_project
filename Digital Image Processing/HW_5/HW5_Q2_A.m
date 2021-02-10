clc
clearvars -except f_fuse;
close all
f=imread('finger.tif');
% figure,imshow(f);

%% generate square image
[M , N] = size(f);
temp = [M N];
ind_min = find(temp == min(temp));
switch ind_min
    case 1
        temp1 = N-M;
        ff = f(1:end , 1+temp1/2 : end - temp1/2);
    case 2
        temp1 = M-N;
        ff = f(1+temp1/2 : end - temp1/2 , 1:end);
end
% figure,imshow(ff)


[M , N] = size(ff);

sigma_x=8;
sigma_y=8;
% T=8;

L_X=6*sigma_x+1;
L_Y=6*sigma_y+1;
LX_2 = fix(L_X/2);
LY_2 = fix(L_Y/2);

tt=0;
T_start = 7; T_step=1; T_end=10;
for T=T_start:T_step:T_end
    p=0; tt= tt+1;
    for theta=0:pi/16:pi-0.01
        p = p+1; 
        for x = -LX_2 : LX_2
            for y = -LY_2: LY_2
                x_theta=x*cos(theta)+y*sin(theta);
                y_theta=-x*sin(theta)+y*cos(theta);
                gb_filter(x+LX_2+1 , y+LY_2+1,p) = (2*pi*sigma_x*sigma_y)^-1 .* exp(-.5*(x_theta.^2/sigma_x^2+y_theta.^2/sigma_y^2)).*cos(2*pi*x_theta/T);
%                 
            end
        end
    end
    for k=1:16
        figure,imshow(gb_filter(:,:,k),[])
        g(:,:,k)=imfilter(ff,gb_filter(:,:,k));
        %     figure, imshow(g(:,:,k)), title(['gabor with k = ',num2str(k)])
    end
    
    %%%%% IMPLmet pixel method to fuse images
    
    [M,N] = size(g(:,:,1));
    W = 5;
    w_2 = floor(W/2);
    std_f = zeros(1,16);
    for m= w_2 +1:W : M - w_2 -1
        for n= w_2 +1:W : N - w_2 -1
            for k=1 : 16
                temp = g(m-w_2 : m+w_2 , n-w_2 : n+w_2,k);
                std_f(1,k) = std2(temp);
            end
            ind_max = min(find( std_f == max(std_f)));
            f_fuse(m-w_2 : m+w_2 , n-w_2 : n+w_2,tt) = g(m-w_2 : m+w_2 , n-w_2 : n+w_2,ind_max);
        end
    end
    
    figure,imshow(f_fuse(:,:,tt),[]), title(['gabor using T =',num2str(T)])
    
end


%%%%%% Fuse METHOD BLOCK
[M,N] = size(f_fuse(:,:,1));
W = 5;
w_2 = floor(W/2);

for T=T_start:T_step:T_end
    std_f = zeros(1,(T_end-T_start)/T_step);
    for m= w_2 +1:W : M - w_2 -1
        for n= w_2 +1:W : N - w_2 -1
            for tt=1 : T_end-T_start +1
                temp = f_fuse(m-w_2 : m+w_2 , n-w_2 : n+w_2,tt);
                std_f(1,tt) = std2(temp);
            end
            ind_max = min(find( std_f == max(std_f)));
            f_fuse_T_A(m-w_2 : m+w_2 , n-w_2 : n+w_2) = f_fuse(m-w_2 : m+w_2 , n-w_2 : n+w_2,ind_max);
        end
    end
end

figure, imshow(f_fuse_T_A,[]),title('final Image fuse with BLOCK method')

%%%%%% Fuse MEthod PIXEL


for T=T_start:T_step:T_end
    std_f = zeros(1,(T_end-T_start)/T_step);
    for m= w_2 +1 : M - w_2 -1
        for n= w_2 +1 : N - w_2 -1
            for tt=1 : T_end-T_start +1
                temp = f_fuse(m-w_2 : m+w_2 , n-w_2 : n+w_2,tt);
                std_f(1,tt) = std2(temp);
            end
            ind_max = min(find( std_f == max(std_f)));
            f_fuse_T_B(m , n) = f_fuse(m , n ,ind_max);
        end
    end
end

figure, imshow(f_fuse_T_B,[]),title('final Image fuse with PIXEL method')