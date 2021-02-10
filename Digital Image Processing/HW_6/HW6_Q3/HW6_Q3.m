clc;
close all;
clear


thr_brown = .13;
thr_white = .2;

% load('brown_mean.mat');
% load('brown_cov.mat');
% 
% load('white_mean.mat');
% load('white_cov.mat');

N0_photo = 10;  %%% Number of phto that we want to use them to extract mean and covariance


%%%%%% Load Images
for i=1:N0_photo
    var_name = sprintf('f%d.jpg',i);
    f{i} = imread(var_name);
end



%%%%% Select Brown points
figure, imshow('select_brown.jpg'),pause(1.5);

for i=1 : N0_photo
    temp =[]; x=[]; y =[];
    while isempty(x)
        figure,imshow(f{i});
        [x,y] = getpts();
        close;
    end
        temp(1,:) = x;
        temp(2,:) = y;
        brown_xy{i} = temp;
end
close;

%%%%% Select White points
figure, imshow('select_white.jpg'),pause(1.5);

for i=1 : N0_photo
temp =[]; x=[]; y =[];
    while isempty(x)
        figure,imshow(f{i});
        [x,y] = getpts();
        close;
    end
    temp(1,:) = x;
    temp(2,:) = y;
    white_xy{i} = temp;
end
close;

%%%%% Calculate mean and covariance matrix
brown_value=[];
white_value=[];

for i=1 : N0_photo
    %%% Extract pixel value of Brown pixels
    temp_brown =  brown_xy{i};
    x_brown = temp_brown(1,:).';
    y_brown = temp_brown(2,:).';
    brown_value = [brown_value; impixel(f{i},x_brown,y_brown)];

    %%% Extract pixel value of White pixels
    temp_white =  white_xy{i};
    x_white = temp_white(1,:).';
    y_white = temp_white(2,:).';
    white_value = [white_value; impixel(f{i},x_white,y_white)];
end

%%% Extract White and Brown mean and covariance

brown_mean = mean(brown_value); %%% mean
white_mean = mean(white_value); %%% mean
for i=1:3
    X(:,:,i) = brown_value(:,i) - brown_mean(1,i);
    Y(:,:,i) = white_value(:,i) - white_value(1,i);
end
brown_cov = [X(:,:,1)' ; X(:,:,2)' ; X(:,:,3)'] * [X(:,:,1) X(:,:,2) X(:,:,3)];
white_cov = [Y(:,:,1)' ; Y(:,:,2)' ; Y(:,:,3)'] * [Y(:,:,1) Y(:,:,2) Y(:,:,3)];




%%%%%%%% we use this code to choose threshode 
% ff = imread('f10.jpg');
% s=1;
% while s==1
%     figure,imshow(ff);
%     [x,y] = getpts();
%     x = impixel(ff,x,y);
%     Mahal_Brown = sqrt((x - brown_mean) * inv(brown_cov) * (x - brown_mean)');
%     Mahal_White = sqrt((x - white_mean) * inv(white_cov) * (x - white_mean)');
%     
%     clc;
%     
%     if Mahal_Brown < Mahal_White && Mahal_Brown < thr_brown
%         disp('chosed pixel is BROWN');
%     elseif Mahal_White < Mahal_Brown && Mahal_White < thr_white
%         disp('chosed pixel is WHITE');
%     else
%         disp('we can''t segment this pixel')
%     end
% end




ff = imread('giraf3_2.jpg');
[M,N] = size(ff(:,:,1));
g_Brown = zeros(M,N);
g_White = zeros(M,N);

for m=1:M
    for n=1:N
        x = impixel(ff,n,m);
        Mahal_Brown = sqrt((x - brown_mean) / (brown_cov) * (x - brown_mean)');
        Mahal_White = sqrt((x - white_mean) / (white_cov) * (x - white_mean)');
        
        if Mahal_Brown < Mahal_White && Mahal_Brown < thr_brown
            g_Brown(m,n) = 1;
        elseif Mahal_White < Mahal_Brown && Mahal_White < thr_white
            g_White(m,n) = 1;
        end
    end
end


% %%%%%%%%%% Black and White

figure,subplot(121),imshow(ff),title('Original Image');
subplot(122),imshow(g_White),title('White Pixel of Giraffe');
suptitle(['Threshold White = ',num2str(thr_white)]) 


figure,subplot(121),imshow(ff),title('Original Image');
subplot(122),imshow(g_Brown),title('Brown Pixel of Giraffe');
suptitle(['Threshold Brown = ',num2str(thr_brown)])

%%%%%%%%% Cloured
[m_b , n_b] = find(g_Brown == 1);
[m_w , n_w] = find(g_White == 1);

Brown_1= zeros(M,N);Brown_2= zeros(M,N);Brown_3= zeros(M,N);
for i=1:length(m_b)
    r = ff(m_b(i),n_b(i),1);
    g = ff(m_b(i),n_b(i),2);
    b = ff(m_b(i),n_b(i),3);
    Brown_1(m_b(i),n_b(i)) = r;
    Brown_2(m_b(i),n_b(i)) = g;
    Brown_3(m_b(i),n_b(i)) = b;
end


White_1= zeros(M,N); White_2= zeros(M,N); White_3= zeros(M,N);
for i=1:length(m_w)
    r = ff(m_w(i),n_w(i),1);
    g = ff(m_w(i),n_w(i),2);
    b = ff(m_w(i),n_w(i),3);
    
    White_1(m_w(i),n_w(i)) = r;
    White_2(m_w(i),n_w(i)) = g;
    White_3(m_w(i),n_w(i)) = b;
end

g_Borwn_Colourd = uint8(cat(3,Brown_1,Brown_2,Brown_3));
g_White_Colourd = uint8(cat(3,White_1,White_2,White_3));


figure,subplot(121),imshow(ff),title('Original Image');
subplot(122),imshow(g_White_Colourd),title('White Pixel of Giraffe');
suptitle(['Threshold White = ',num2str(thr_white)]) 


figure,subplot(121),imshow(ff),title('Original Image');
subplot(122),imshow(g_Borwn_Colourd),title('Brown Pixel of Giraffe');
suptitle(['Threshold Brown = ',num2str(thr_brown)])