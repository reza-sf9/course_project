clc;
clear;
close all;

% load('f_skeleton.mat');
f=imread('1.png');
f= im2double(f);
load('BW.mat','BW');    %%% Load Strels
[M,N] = size(f);

figure,imshow(f),title('Skeleton Image');

non_zero_skeletion = length(find(f)); %%%  # of non_zero elements of skeletion (the stopping criteria is based on this)

counter = 0; %%%  # of itteration of applying strels on Image
g=f;
exit_1 =0; %%% stopping cirteria of Main Loop
while exit_1 == 0
    counter = counter + 1;
    for i=1:8
        BBW = BW(:,:,i); %%% using strek No i
        w = length(BBW);
        w_2 = floor(w/2);
        exit_2 = 0;
        
        while exit_2 == 0
            
            ff = zeros(M,N);
            
            ind_neg = find(BBW == -1); l_ind_neg = length(ind_neg); %%% indices that must be zero
            ind_pos = find(BBW ==  1); l_ind_pos = length(ind_pos); %%% indices that must be one
            
            for m = 1+w_2 : M-w_2
                for n = 1+w_2 : N-w_2
                    
                    temp = f(m-1:m+1,n-1:n+1); %%% window that we want compare it with strel
                    temp_neg = find(temp == 0);  l_temp_neg = length(temp_neg); %%% indices that r 1 of window (center pixel is m,n)
                    temp_pos = find(temp ==  1);  l_temp_pos = length(temp_pos);%%% indices that r 0 of window (center pixel is m,n)
                    
                    %%%% Assessing weather 1 indices are exit or not
                    c1=0;
                    for j =1:l_ind_neg
                        c1 = c1 + length(find(ind_neg(j) == temp_neg));
                    end
                    
                    c2=0;
                    if c1 == l_ind_neg
                        %%%% Assessing weather 1 indices are exit or not
                        for j=1:l_ind_pos
                            c2 = c2 + length(find(ind_pos(j) == temp_pos));
                        end
                        %%%% if 0 & 1 indices are exist, Hit-Miss ocures
                        if c2 == l_ind_pos
                            ff(m,n) =1;
                        end
                    end
                    
                end
            end
            
            f = f - ff;
            aa_2 = find(ff == 1);
            if isempty(aa_2) %%%% if this condition is True , we have No change with this strel and must change or strel
                exit_2 =1;
            end
        end
    end
    
    f = f-ff;
    l_f = length(find(f==1));
    l_g = length(find(g==1));
    if  l_f == l_g %%%% if this condition is True , we have No change and we must exit
        exit_1 =1;
    end
    g=f;
    tit = sprintf('Thining of Skelton Image -- No of itteration = %d \n No of changes pixel this itteration = %d',counter,l_g-l_f);
    figure,imshow(f),title(tit)
end

