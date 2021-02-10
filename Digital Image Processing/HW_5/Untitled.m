clc
clear 
close all

f=imread('photo561829177866693988.jpg');
f = rgb2gray(f);
f = im2double(f);
[M,N] = size(f);

% g = zeros(M,N);
% 
% thr = graythresh(f)+.1;
% ind_1 = find(f>= thr);
% g(ind_1) = 1;
% figure, imshow(g)


W=100;
ff = zeros(M,N);
gg = zeros(W,W);
for i=1:W:M-W
    for j=1:W:N-W
        gg = zeros(W,W);
       temp = f(i:i+W-1,j:j+W-1);
       thr = graythresh(f)+.1;
       ind_11 = find(temp>= thr);
       gg(ind_11) = 1;
       figure,imshow(gg),title(['i= ',num2str(i),'  j= ',num2str(j)]);
       ff(i:i+W-1,j:j+W-1) = gg;
    end
end
figure,imshow(ff);