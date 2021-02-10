clc
close all
clear

r_max = 200*sqrt(2); r_min = 0;
T_max = 40;          T_min = 10;

p=0; 
for m=-200:200
    p=p+1; q=0;
    for n=-200:200
        q=q+1;
        r= sqrt(m^2 + n^2);
        T = ((T_min -T_max)/(r_max - r_min))*r + T_max;
        f(p,q) = sin(2*pi*r/T);
    end
end

figure 
imshow(f,[])
title('Original Picture');
