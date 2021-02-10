clc;
clear;
close all;

%%%%%%%%%
f = imread('baby2.jpg');
f= rgb2gray(f);
f = im2double(f);
figure, imshow(f) , title('Oroginal Image');
[M,N] = size(f);
thr = .5;
%%%%%%%    error difusion
h_error = [0 0 7/16;3/16 5/16 1/16];

f_error = f;
for i=1 : M -2
    for j=2: N -2

        temp_1 = f_error(i:i+1 , j-1:j+1);
        temp_2 = f_error(i,j);

        if temp_2<= thr
            f_error(i,j) = 0;
        else
            f_error(i,j) = 1;
        end
        e = f_error(i,j) - temp_2;
        temp_3 = e.*h_error;
        f_error(i:i+1 , j-1:j+1) = f_error(i:i+1 , j-1:j+1) - temp_3;
    end
end

figure, imshow(f_error) , title('Error Difusion')

%%%%%%%%%%  Dot Fussion Optimized Method

C_Mat = load('C_Matrix.txt');

f_dot= zeros(8*ceil(M/8),8*ceil(N/8));
f_dot(1:M,1:N) = f ;
[M1 , N1] = size(f_dot);
temp_1 = zeros(10,10);
C_Mat_1 = zeros(10,10);
C_Mat_1(2:9,2:9) = C_Mat;
h_dot = zeros(3,3);
p=0;
for i=1:8: M1
    for j=1:8: N1
        p=p+1;
        temp_1(2:9,2:9)=f_dot(i:i+7,j:j+7);
        for k=1:64
            [row , col] = find(C_Mat_1 ==k);
            
            if temp_1(row,col) > thr
                e = 1 - temp_1(row,col);
                temp_1(row,col) = 1;
            else
                e = 0 - temp_1(row,col);
                temp_1(row,col) = 0;
            end
            
            temp_2 = temp_1(row-1:row+1 , col-1: col+1);
            temp_3 = C_Mat_1(row-1:row+1 , col-1: col+1);
            
            index_1 = find(temp_3 > k);           %%% indeices that greater than k
            if isempty(index_1) == 0
                index_1_e = index_1(mod(index_1,2)==0);
                index_1_o = index_1(mod(index_1,2)==1);
                h_dot(index_1_e) = 2;
                h_dot(index_1_o) = 1;
                
                index_2 = find(temp_3 <= k);   %%% indeices that greater than k
                h_dot(index_2) = 0;
                h_dot = (1/sum(h_dot(:))).*h_dot;
                temp_4 = temp_2 - e.*h_dot;
                temp_1(row-1:row+1 , col-1: col+1) = temp_4;
            end
            
        end
        f_dot(i:i+7,j:j+7) = temp_1(2:9,2:9);
        
        
    end
end
f_dot = f_dot(1:M,1:N);
figure, imshow(f_dot),title('dot difussion new Method')

%%%%%%%%%%  Dot Fussion Knuth Method

C_Mat = load('C_Matrix.txt');
Knuth_Mat = load('Knuth_Mat.txt');

f_dot= zeros(8*ceil(M/8),8*ceil(N/8));
f_dot(1:M,1:N) = f ;
[M1 , N1] = size(f_dot);
temp_1 = zeros(10,10);
Knuth_Mat_1 = zeros(10,10);
Knuth_Mat_1(2:9,2:9) = Knuth_Mat;
h_dot = zeros(3,3);
p=0;
for i=1:8: M1
    for j=1:8: N1
        p=p+1;
        temp_1(2:9,2:9)=f_dot(i:i+7,j:j+7);
        for k=1:64
            [row , col] = find(Knuth_Mat_1 ==k);
            
            if temp_1(row,col) > thr
                e = 1 - temp_1(row,col);
                temp_1(row,col) = 1;
            else
                e = 0 - temp_1(row,col);
                temp_1(row,col) = 0;
            end
            
            temp_2 = temp_1(row-1:row+1 , col-1: col+1);
            temp_3 = Knuth_Mat_1(row-1:row+1 , col-1: col+1);
            
            index_1 = find(temp_3 > k);           %%% indeices that greater than k
            if isempty(index_1) == 0
                index_1_e = index_1(mod(index_1,2)==0);
                index_1_o = index_1(mod(index_1,2)==1);
                h_dot(index_1_e) = 2;
                h_dot(index_1_o) = 1;
                
                index_2 = find(temp_3 <= k);   %%% indeices that greater than k
                h_dot(index_2) = 0;
                h_dot = (1/sum(h_dot(:))).*h_dot;
                temp_4 = temp_2 - e.*h_dot;
                temp_1(row-1:row+1 , col-1: col+1) = temp_4;
            end
            
        end
        f_dot(i:i+7,j:j+7) = temp_1(2:9,2:9);
        
        
    end
end
f_dot = f_dot(1:M,1:N);
figure, imshow(f_dot),title('dot difussion Knuth Method')