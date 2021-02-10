clc;
clear;
close all;

I = imread('glass2.jpg');
I = rgb2gray(I);

figure,imshow(I), title('Original Image');

%%%% Add salt and pepper noise
f = imnoise(I,'salt & pepper',0.4);
figure, imshow(f),title('Image with Salt and Pepper noise');

[M,N] = size(f);



temp_I = f;
W=9; %%% max length of Window
w=3; %%% initial value of w
l_w = floor(w/2);
l_W = floor(W/2);

p=0;
for m=l_W+1 : M - l_W
    p= p+1; q=0;
    for n=l_W+1 : N - l_W
        l_w = floor(w/2);
        q=q+1;
        ch =1;
        while ch == 1;
            ff = temp_I(m,n);
            
            temp = temp_I(m-l_w : m+l_w , n-l_w : n+l_w);
            sort_t = sort(temp(:));
            
            med_t = median(sort_t); 
            min_t = min(sort_t); 
            max_t = max(sort_t);
            
            if min_t< med_t && max_t> med_t
                if min_t< ff && max_t> ff
                    g(p,q) = temp_I(m,n);
                    ch=0;
                else
                    g(p,q) = med_t;
                    ch=0;
                end
            elseif w<W
                w= w+2;
                l_w = floor(w/2);
                ch=1;
            else
                g(p,q) = ff;
                ch=0;
            end
        end
        
    end
end

gcf = figure;
imshow(g),title('denoising with Adaptive Median Filter method')


