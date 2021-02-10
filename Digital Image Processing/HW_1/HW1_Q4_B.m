clc;
clear ; 
close all;

%%%%%%%%%%%%%% 3
teta_max = 90;       teta_min = 0;
T_max = 30;          T_min = 8;

for m=0 : 150
    for n= 0 :150
        r= sqrt(m^2 + n^2);
        teta = atand(m/n);
%         teta(isnan(teta))=0;
        
        T = ((T_min -T_max)/(teta_max - teta_min))*teta + T_max;
        f(m+1,n+1) = sin(2*pi*r/T);
    end
end

figure 
imshow(f,[])
title('cute teammate');

%%%%%%%%%%%%%%%%%%% 4
% f=f>=0;
%%%% nearest neighbohood interpolation
S = 3.5;
step = 1/S;
p = 0; 
for m=0 :step : 150
    p = p+1; q = 0;
    for n= 0 :step :150
        q = q+1;
         l= round(m);
         k=round(n);

        g(p,q) = f(l+1,k+1);
    end
end

figure 
imshow(g,[])
str_tit = sprintf('nearest neighbor interpolation S=%d ',S);
title(str_tit);

%%%%% linear imterpolation 
S = 2.5;
step = 1/S ;

p = 0; 
for m=0 :step : 149
    p = p+1; q = 0;
    for n= 0 :step :149
        q = q+1;
        
         l= floor(m)+1;
         k= floor(n)+1;
         
         a = f(l , k) ; 
         b = f(l , k+1);
         c = f(l+1 , k);
         d = f(l+1 , k+1);
         
         A = step*c + (1-step)*a;
         B = step*d + (1-step)*b;
 
        h(p,q) = step*B + (1-step)*A;
    end
end

figure 
imshow(h,[])
str_tit = sprintf('linear interpolation S=%d ',S);
title(str_tit);



