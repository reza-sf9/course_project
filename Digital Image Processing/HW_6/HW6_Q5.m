clc;
clearvars -except f_local_b;
close all;

f = imread('finger.tif');
ff = im2double(f);




% figure,imshow(f),title('Oroginal Image');

[M,N] = size(ff);

W = 7; w_2 = floor(W/2);
f_local_b = zeros(M,N);
for m=(w_2+1) :5: M - (w_2+1)
    for n= (w_2+1) :5: N - (w_2+1)
        temp = ff(m-w_2 : m+w_2, n-w_2 : n+w_2);
        thr = graythresh(temp);
        if thr ~=0
           k=0; 
        end
        temp_2 = zeros(W,W);
        for mm=1:W
            for nn=1:W
               if   temp(mm,nn) > thr
                    temp_2(mm,nn) = 1;     
               end
            end
        end
        f_local_b(m-w_2 : m+w_2, n-w_2 : n+w_2) = temp_2;
    end
end

h=figure,imshow(f_local_b),title(['Local Binrization with OTSU Method -- Win Length = ' , num2str(W)]);
tit = sprintf('Binrization_W_%d.jpg',W);
saveas(h,tit);


r1 =2;
strl = strel('disk',r1,0);
f1 = imclose(f_local_b , strl);
h1 = figure,imshow(f1),title(['Image after imclose with a disk with radius = ',num2str(r1)]);
tit_1 = sprintf('W_%d_imclose_Disk_%d.jpg',W,r1);
saveas(h1,tit_1);


r2 = 2;
strl = strel('square',r2);
f2 = imopen(f1 , strl);
h2 = figure,imshow(f2),title(['Image after imopen with a square with side = ',num2str(r2)]);
tit_1 = sprintf('W_%d_imopen_Square_%d.jpg',W,r2);
saveas(h1,tit_1);






 %%%%% Binrization with 1 threshold
thr = graythresh(ff).*255;
ff = ff.*255;


f_b = zeros(M,N);
for i=0:255
    ind = find(ff == i);
    if i > thr
       f_b(ind) = 1;
    end
end
figure,imshow(f_b),title('Binarization Image with Global THreshold');