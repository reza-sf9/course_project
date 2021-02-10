clc;
clear;
close all;

syms x;
figure,ezplot(log(exp(x) + exp(-x)), [-10, 10])

figure,ezplot(-log(x) + x, [-1, 10])


I = imread('glass2.jpg');
I = rgb2gray(I);

% figure,imshow(I), title('Original Image');

%%%% Add salt and pepper noise
f = imnoise(I,'salt & pepper',0.4);
figure, imshow(f),title('Image with Salt and Pepper noise');

[M,N] = size(f);


for W=3:-2:3                 %% window length
    temp_I = f;
    l_w = floor(W/2);
    for itt=1:8              %% # of itteration
        p=0; g=[];
        for m=l_w+1 : M - l_w
            p= p+1; q=0;
            for n=l_w+1 : N - l_w
                q=q+1;
                temp = temp_I(m-l_w : m+l_w , n-l_w : n+l_w);
                sort_t = sort(temp(:));
                med = median(sort_t);
                
                g(p,q) = med;
            end
        end
        temp_I(l_w+1 : M - l_w,l_w+1 : N - l_w) = g;
        gcf = figure;
        imshow(g,[]),title(['denoising with lenght of window =',num2str(W),' and # of ittereation ',num2str(itt)])
        tit_pic = sprintf('denoise_W_%d___ITT_%d.jpg',W,itt);
        saveas(gcf,tit_pic)
    end
end

% figure, imshow(g),title(['denoising with lenght of window =',num2str(W),' and # of ittereation ',num2str(itt)])