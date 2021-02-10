clc;
clear;
close all;
%%%% Generate chirp picture
r_max = 200*sqrt(2); r_min = 0;
T_max = 50;          T_min = 10;

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
f_Binary = f>=0;   %%%% convert photo 2 binary  

 figure  
 imshow(f,[])
 title('Chirp Picture');
% %%% Binary show
 figure  
 imshow(f_Binary,[])
 title('Binary show of Chirp Picture');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   downsample photo 
S =2.3; %%% downsample rate
[M,N] = size(f);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% nearest neighbohood interpolation
%%%%%%%%%% DOWNSAMPLE
p = 0; 
for m=1 :S : M
    p = p+1; q = 0;
    for n= 1 :S :N
        q = q+1;
         l= round(m);
         k=round(n);
        f_Nearest_D(p,q) = f(l,k);
    end
end
f_Nearest_D_B = f_Nearest_D>=0;
 figure 
 imshow(f_Nearest_D_B,[])
 str_tit = sprintf('nearest neighbor decrease with S=%.1f ',S);
 title(str_tit);
%%%%%%%% UPSAMPLE
[M,N] = size(f_Nearest_D);
p = 0; 
for m=1 : 1/S : M
    p = p+1; q = 0;
    for n= 1 :1/S :N
        q = q+1;
         l= round(m);
         k=round(n);
        f_Nearest_U(p,q) = f_Nearest_D(l,k);
    end
end
f_Nearest_U_Binary=f_Nearest_U>=0;    %%% convert photo to a binary photo 

figure 
imshow(f_Nearest_U,[])
str_tit = sprintf('Nearest Neighbor increase with S=%.1f ',S);
title(str_tit);
%%%%% Binary show 
figure 
imshow(f_Nearest_U_Binary,[])
str_tit = sprintf('Binaryb show of Nearest Neighbor interpolation with S=%.1f ',S)
title(str_tit)

[M,N] = size(f_Nearest_U);
error_Neighbor = norm(f(1:M,1:N)-f_Nearest_U)/norm(f)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% linear interpolation
%%%%%%%% DOWNSAMPLE
[M,N] = size(f);
p = 0; 
for m=1 :S : M
    p = p+1; q = 0;
    for n= 1 :S :N
        q = q+1;
        
        step_m = m-floor(m);
        step_n = n-floor(n);
        
        
         l= floor(m);
         k= floor(n);
         
         a = f(l , k) ; 
         b = f(l , k+1);
         c = f(l+1 , k);
         d = f(l+1 , k+1);
         
         A = step_m*c + (1-step_m)*a;
         B = step_m*d + (1-step_m)*b;
 
        f_Linear_D(p,q) = step_n*B + (1-step_n)*A;
    end
end
f_Linear_D_B = f_Linear_D>=0;
 figure 
 imshow(f_Linear_D_B,[])
 str_tit = sprintf('linear decrease with S=%.1f ',S);
 title(str_tit);

%%%%%%% UPSAMPLE
[M,N] = size(f_Linear_D);
p = 0; 
for m=1 : 1/S : M-1
    p = p+1; q = 0;
    for n= 1 : 1/S :N-1
        q=q+1;
        
        step_m = m-floor(m);
        step_n = n-floor(n);
        
         l= floor(m);
         k= floor(n);
         
         a = f_Linear_D(l , k) ; 
         b = f_Linear_D(l , k+1);
         c = f_Linear_D(l+1 , k);
         d = f_Linear_D(l+1 , k+1);
         
         A = step_m*c + (1-step_m)*a;
         B = step_m*d + (1-step_m)*b;
 
        f_Linear_U(p,q) = step_n*B + (1-step_n)*A;
    end
end
f_Linear_U_Binary=f_Linear_U>=0;    %%% convert photo to a binary photo 

figure 
imshow(f_Linear_U,[])
str_tit = sprintf('linear  increase  with S=%.1f ',S);
title(str_tit);
%%%%% Binary show
figure;
imshow(f_Linear_U_Binary,[])
str_tit = sprintf('Binary show of linear interpolation with S=%.1f ',S);
title(str_tit);

[M,N] = size(f_Linear_U);
error_Linear = norm(f(1:M,1:N)-f_Linear_U)/norm(f)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% spline Interpolation
A =[-0.5 1.5 -1.5 0.5; 1 -2.5 2 -0.5;-0.5 0 0.5 0;0 1 0 0];
p_c = (3:-1:0).';                              %%% power coef
%%%%%%% DOWN SAMPLE
[M,N] = size(f);
f_ZP = zeros(M+3,N+3);                         %%% zero-pad a row and column at the first and 2 row and coloumn at the end
f_ZP(2:end-2,2:end-2) = f; 
[M,N] = size(f_ZP);
p=0;
for m=2:S:M-2
    p=p+1; q=0;
    for n=2:S:N-2
        q=q+1;
        f_m = fix(m);                           %%% fixed value of m 
        d_m = m - f_m;                          %%% delta_m
        f_n = fix(n);                           %%% fixed value of n 
        d_n = n - f_n;                          %%% delta_n
        
        IP_Data = f_ZP(f_m-1:f_m+2,f_n-1:f_n+2);   %%% extract 16  points in order to interpolate
        
        %%%% interpolatin for 1st row
        IP_A = IP_Data(1,:).';                  %%% interpolate with 4 points
        Coef_A = A*IP_A;                        %%% calculate coeficients
        IP_ABCD(1,1) = Coef_A.' * (d_n.^ p_c);  %%% calulate interpolated point
        
        %%%% interpolatin for 2nd row
        IP_B = IP_Data(2,:).';
        Coef_B = A*IP_B;                        
        IP_ABCD(2,1) = Coef_B.' * (d_n.^ p_c);  
        
        %%%% interpolatin for 3rd row
        IP_C = IP_Data(3,:).';
        Coef_C = A*IP_C;                        
        IP_ABCD(3,1) = Coef_C.' * (d_n.^ p_c);  
        
        %%%% interpolatin for 4th row
        IP_D = IP_Data(4,:).';
        Coef_D = A*IP_D;                        
        IP_ABCD(4,1) = Coef_D.' * (d_n.^ p_c);  
        
        %%%% interpolatin of 4 above row
        Coef_Tot = A*IP_ABCD;                        
        f_A_Spline_D(p,q) = Coef_Tot.' * (d_m.^ p_c);  %%% f down sample with spline approach
    end
end
f_A_Spline_D_B = f_A_Spline_D>=0;
 figure 
 imshow(f_A_Spline_D_B,[])
 str_tit = sprintf('cubic decrease with S=%.1f ',S);
 title(str_tit);

%%%%%% UPSAMPLE
[M,N] = size(f_A_Spline_D);
f_A_Spline_D_ZP = zeros(M+3,N+3);                         %%% zero-pad a row and column at the first and 2 row and coloumn at the end
f_A_Spline_D_ZP(2:end-2,2:end-2) = f_A_Spline_D; 
[M,N] = size(f_A_Spline_D_ZP);

p=0;
for m=2:1/S:M-2
    p=p+1; q=0;
    for n=2:1/S:N-2
        q=q+1;
        f_m = fix(m);                           %%% fixed value of m 
        d_m = m - f_m;                          %%% delta_m
        f_n = fix(n);                           %%% fixed value of n 
        d_n = n - f_n;                          %%% delta_n
        
        IP_Data = f_A_Spline_D_ZP(f_m-1:f_m+2,f_n-1:f_n+2);   %%% extract 16  points in order to interpolate
        
        %%%% interpolatin for 1st row
        IP_A = IP_Data(1,:).';                  %%% interpolate with 4 points
        Coef_A = A*IP_A;                        %%% calculate coeficients
        IP_ABCD(1,1) = Coef_A.' * (d_n.^ p_c);  %%% calulate interpolated point
        
        %%%% interpolatin for 2nd row
        IP_B = IP_Data(2,:).';
        Coef_B = A*IP_B;                        
        IP_ABCD(2,1) = Coef_B.' * (d_n.^ p_c);  
        
        %%%% interpolatin for 3rd row
        IP_C = IP_Data(3,:).';
        Coef_C = A*IP_C;                        
        IP_ABCD(3,1) = Coef_C.' * (d_n.^ p_c);  
        
        %%%% interpolatin for 4th row
        IP_D = IP_Data(4,:).';
        Coef_D = A*IP_D;                        
        IP_ABCD(4,1) = Coef_D.' * (d_n.^ p_c);  
        
        %%%% interpolatin of 4 above row
        Coef_Tot = A*IP_ABCD;                        
        f_A_Spline_U(p,q) = Coef_Tot.' * (d_m.^ p_c);  %%% f up sample with spline approach
    end
end
f_A_Spline_U_Binary=f_A_Spline_U>=0;    %%% convert photo to a binary photo 

figure 
imshow(f_A_Spline_U,[])
str_tit = sprintf('cubic increase with S=%.1f ',S);
title(str_tit);
%%% Binary show
figure
imshow(f_A_Spline_U_Binary,[])
str_tit = sprintf('Binary show of cubic interpolation with S=%.1f ',S)
title(str_tit);

[M,N] = size(f_A_Spline_U);
error_Cubic = norm(f(1:M,1:N)-f_A_Spline_U)/norm(f)

%%%% comparison 1D f , linear  spline and Nearst
U_point = 390;
a=f(1,1:U_point);
b=f_A_Spline_U(1,1:U_point);
c=f_Linear_U(1,1:U_point);
d=f_Nearest_U(1,1:U_point);
figure
plot(d,'b'),hold on;         %% Nearest
plot(c,'k:'),hold on;          %% Linear
plot(b,'g'),hold on;           %% spline
plot(a,'r');  %% original f
