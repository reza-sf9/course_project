clc
clear
close all;

g=imread('bloodcel_95.jpg');
f=rgb2gray(g);
f = double(f);
[M,N] = size(f);
ff = ones(M,N);
ind_x = find(f<100);
ff(ind_x) = 0;    %binerization of image
figure, imshow(ff);
title('binerization of cells image with threshold')

% making strel and using it for hit and miss
r=13;
W = 2*r+1;
strl = zeros(W , W);
for m=-r:r
    for n=-r:r
        d = sqrt(m.^2 + n.^2);
        if ceil(d) <= r
            strl(m+r+1,n+r+1) = -1;
        end
    end
end
figure,imshow(strl,[]),title(' disk strel with r=13');
BBW = strl;     %%% using strel


g = ones(M+W , N+W);
g(1+r:M+r , 1+r:N+r) = ff;
[MM,NN] = size(g);



gr = bwhitmiss(ff,BBW);
gr = double(gr);
figure,imshow(gr)
title('hit and miss binary image with disk strel');



p=0;
cond=0;
grr = zeros(M,N);
[ind_x,ind_y] = find(gr == 1);
thr_r = 13;
%this loop is for choosing one 1 for every cell
while cond==0
    p=p+1;
    
    r = sqrt( (ind_x(1:end)-ind_x(1)).^2 + (ind_y(1:end)-ind_y(1)).^2 );
    h = find(r < thr_r);
    ind_n_x = round(mean(ind_x(h)));
    ind_n_y = round(mean(ind_y(h)));
    
    grr(ind_n_x,ind_n_y) =1;
    for i=1:length(h)
        gr(ind_x(h(i)),ind_y(h(i))) = 0;
    end
    
    [ind_x,ind_y] = find(gr == 1);
    if isempty(ind_x)
        cond =1;
    end
    
end

figure,imshow(grr)
[ind_xx , ind_yy] = find(grr==1);

number_cell = length(ind_xx)     %calculating the numer of ones as cell numbers


