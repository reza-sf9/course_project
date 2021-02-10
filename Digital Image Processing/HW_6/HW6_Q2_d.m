clc;
clear;
close all;

load('f_WLennght_23_Thr_.105.mat');
[M,N] = size(g);
R=15;
f = zeros(M+2*R , N+2*R);
f(R+1:end-R, R+1:end-R)=g;
figure,imshow(f);
title('the edges of cells');


thr = .5;
temp_R = R.* ones(2*R+1,2*R+1);


fff = zeros(M+2*R,N+2*R);% 
[M_1,N_1] = find(f); %%% find indice that are equal to 1

ff= zeros(M+2*R , N+2*R);
for i=1:length(M_1)
    p=0;
   for m=M_1(i) - R : M_1(i) + R
       q=0;p=p+1;
       for n=N_1(i) - R : N_1(i) + R
           q=q+1;
           r(p,q) = sqrt((m-M_1(i))^2 + (n-N_1(i))^2);
       end
   end
  dif_temp = abs(r - temp_R);
  ind_1 = find(dif_temp < thr); %% these inice must be 1
  temp_f = zeros(2*R+1,2*R+1);
  temp_f(ind_1) = 1;
  temp_ff = zeros(M+2*R , N+2*R);
  jj=M_1(i) - R;
  jjj=M_1(i) + R ;
  jjjj = M_1(i);
  temp_ff(M_1(i) - R:M_1(i) + R  , N_1(i) - R:N_1(i) + R ) = temp_f;
  ff = ff + temp_ff;
end

ind = find(ff>15);
fff(ind) = 1;
[NN,MM] = meshgrid(1:N+2*R,1:M+2*R);
figure,mesh(MM,NN,fff),xlabel('X') , ylabel('Y');
title('showing the places of cells with the center of them');

p=0;
cond=0;
g = zeros(M+2*R,N+2*R);
[ind_x,ind_y] = find(fff == 1);

while cond==0
    p=p+1;
    
    r = sqrt( (ind_x(1:end)-ind_x(1)).^2 + (ind_y(1:end)-ind_y(1)).^2 );
    h = find(r < 4.5*R);
    ind_n_x = round(mean(ind_x(h)));
    ind_n_y = round(mean(ind_y(h)));
    
    g(ind_n_x,ind_n_y) =1;
    for i=1:length(h)
        fff(ind_x(h(i)),ind_y(h(i))) = 0;
    end
    
    [ind_x,ind_y] = find(fff == 1);
    if isempty(ind_x)
        cond =1;
    end
    
end

figure,imshow(g)
[ind_xx , ind_yy] = find(g==1);

%%%% Create Circles
gg=g;
for i=1:length(ind_xx)
   ind_x = ind_xx(i);
   ind_y = ind_yy(i);
   p=0; r=[];
   for m=ind_x - R:ind_x + R 
       q=0;p=p+1;
       for n=ind_y - R:ind_y + R 
           q=q+1;
           r(p,q) = sqrt((m-ind_x)^2 + (n-ind_y)^2); 
       end
   end
  dif_temp = abs(r - temp_R);
%   ind_1 = find(dif_temp <= thr); %% these inice must be 1
  ind_1 = find(r <= R); %% these inice must be 1
  temp_f = zeros(2*R+1,2*R+1);
  temp_f(ind_1) = 1;
  gg(ind_x - R:ind_x + R  , ind_y - R:ind_y + R ) = temp_f;
end
gg = gg(R+1:end-R,R+1:end-R);
figure,imshow(gg),title('the circles with the R=15');

[M_ ,N_] = size(gg);
ggg = zeros(M_ , N_);
ind_0 = find(gg ==0);
ggg(ind_0) = 1;
figure,imshow(ggg),title('ggg');

f = imread('bloodcel_95.jpg');
f = rgb2gray(f);
f = im2double(f);
figure,imshow(f), title('Original Image')
