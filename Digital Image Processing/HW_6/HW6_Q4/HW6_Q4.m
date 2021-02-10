clc;
clearvars -except xx yy f_n;
% clear;
close all;

f = 25.*ones(200,200);
%%%%%% Initial values
RR=[30 20 15];
xx_center = [70 70 140];
yy_center = [70 140 100];
incriment = [10 5 25];
%%%% Create Image

for i=1:3
    R = RR(i);
    x_center = xx_center(i);
    y_center = yy_center(i);
    
    %%% Create circles
    if i==1 || i==2
        temp_R = zeros(2*R+1,2*R+1);
        for m=-R:R
            for n=-R:R
                if sqrt(m^2 + n^2) <= R
                    temp_R(m+R+1,n+R+1) = incriment(i);
                end
            end
        end
        ff = zeros(200,200);
        ff(x_center-R:x_center+R,y_center-R:y_center+R) = temp_R;
        f = ff + f;
        
        %%%%%% Create Square
    else
        temp_R = incriment(i) .* ones(2*R+1,2*R+1);
        ff = zeros(200,200);
        ff(x_center-R:x_center+R,y_center-R:y_center+R) = temp_R;
        f = ff + f;
    end
end
figure,imshow(f,[]),title('Original Image');

%%%% Generate Noise
n=2*randn(200,200);
f_n=f+n;
figure,imshow(f_n,[]),title('f + Noise');

%%%%% Create Binary Image
f_Binary = zeros(200,200);
ind_1 = find(f==25);
ind_2 = find(f>25);
f_Binary(ind_1) = 0;
f_Binary(ind_2) = 1;
 figure,imshow(f_Binary,[]),title('Binary f');
[M,N] = size(f);

%%%%%%%%%%%%%%%%

%%%%%%%%%%%% A
% [yy,xx] = get_points(f_n);

alpha = 5; beta=25; gamma = 250; win_length = 7;
Kass_Method(f_n , yy , xx , alpha , beta , gamma , win_length);

%%%%%%%%%%%% B
Level_Set_Method(f)

%%%%%%%%%%% C
Morphology_Method(f_Binary)
% 
% %%%%%%%%%%%% D
f=imrotate(f_Binary,240,'bilinear','crop');
squre_Detection(f)



