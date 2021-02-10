clc;
clear;
close all;

load('exersize4_2.mat');
original_code = round(z);

f=imread('lena.jpg');

%%%%% first we must put our code in a vector
[M,N] = size(f);
code_length = (M*N/64);
code = zeros(1,code_length); %%% we must do zeropad , at the end we will remove additional zeros
temp =original_code.';
code(1,1:20*97) = temp(:);


%%%%
f_new = zeros(M,N);
p=0; counter_code =0;

%%%%%% selecting 2 points in order to coding process
pix_1_x = 4; pix_1_y = 1;
pix_2_x = 2; pix_2_y = 3;

for m=1:8:M
    p=p+1; q=0;
    for n=1:8:N
        q=q+1; counter_code = counter_code+1;
 
        DCT_temp = dct2(f(m:m+7,n:n+7));
        
        pixel_code = code(counter_code);
        
        switch pixel_code
            case 0
                if DCT_temp(pix_1_x,pix_1_y) >= DCT_temp(pix_2_x,pix_2_y)
                    temp = DCT_temp(pix_1_x,pix_1_y);
                    DCT_temp(pix_1_x,pix_1_y) = DCT_temp(pix_2_x,pix_2_y);
                    DCT_temp(pix_2_x,pix_2_y) = temp;
                end
            case 1
                if DCT_temp(pix_1_x,pix_1_y) <= DCT_temp(pix_2_x,pix_2_y)
                    temp = DCT_temp(pix_1_x,pix_1_y);
                    DCT_temp(pix_1_x,pix_1_y) = DCT_temp(pix_2_x,pix_2_y);
                    DCT_temp(pix_2_x,pix_2_y) = temp;
                end
        end
        
       f_new(m:m+7,n:n+7) = idct2(DCT_temp);
    end
end


figure,subplot(121),imshow(f),title('Original Image');
subplot(122),imshow(f_new,[]),title(['Image with DCT Steganography ,Using pixel(' num2str(pix_1_x),','...
    ,num2str(pix_1_y),') & pixel(',num2str(pix_2_x),',',num2str(pix_2_y),')']);

%%%%% decoding

[M,N] = size(f_new);
index_1 = zeros(1,code_length); index_2 = zeros(1,code_length);
q=0;
for m=1:8:M
    for n=1:8:N
        q=q+1;
   temp = f_new(m:m+7,n:n+7); 
   DCT_temp = dct2(temp);
   index_1(1,q) = DCT_temp(pix_1_x,pix_1_y);        
   index_2(1,q) = DCT_temp(pix_2_x,pix_2_y);     
    end
end

temp_code = zeros(1,code_length);
index_1 = find(index_1 - index_2 > 0); %%% find index that the code is 1
temp_code(index_1) = 1;

vector_code = temp_code(1,1:20*97);         %%% delete the zeropad that we added at the first 
retrieve_code = reshape(vector_code ,[97 20]);

% gg = original_code - retrieve_code;

figure,subplot(211),imshow(original_code),title('Oroginal Message');
subplot(212),imshow(retrieve_code.',[]),title(['Retrieve Message, using Pixel(' num2str(pix_1_x),',',num2str(pix_1_y)...
    ') & pixel(',num2str(pix_2_x),',',num2str(pix_2_y),')'])
