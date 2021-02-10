clc;
clear
close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%% without chirp
p=0; 
for m=-200:200
    p=p+1; q=0;
    for n=-200:200
        q=q+1;
        r= sqrt(m^2 + n^2);
        T=10;
        f(p,q) = sin(2*pi*r/T);
    end
end

figure
imshow(f,[])
title('Using a Single Tone Sine with T=10');

ff= fft2(f);
step=2*pi/length(ff);
[U,V]=meshgrid(0:step:2*pi-step , 0:step:2*pi-step);
figure
mesh(U,V ,abs(ff))
title('Fourier Transform of Single Tone Sine')
xlabel('U(radian)'), ylabel('V(radian)')

%%%%%%%% down sample
S=3;    %% scale of downsampling

p=0; 
for m=-200:S:200
    p=p+1; q=0;
    for n=-200:S:200
        q=q+1;
        g(p,q) = f(m+201,n+201);
    end
end

figure 
imshow(g,[])
str_tit = sprintf('Sigle Tone picture after downsampling with S=%d (without chirp)',S);
% title(str_tit);

gg= fft2(g);
step=2*pi/length(gg);
[U,V]=meshgrid(0:step:2*pi-step , 0:step:2*pi-step);
figure
mesh(U,V ,abs(gg))
str_tit = sprintf(' Fourier Transfporm of Sigle Tone picture after downsampling with S=%d (without chirp)',S);
title(str_tit);
xlabel('U(radian)'), ylabel('V(radian)')




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% chirp
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
title('Original Picture(chirp)');

ff= fft2(f);
step=2*pi/length(ff);
[U,V]=meshgrid(0:step:2*pi-step , 0:step:2*pi-step);
figure
mesh(U,V ,abs(ff))
title('Fourier Transform of Original Picture(chirp)')
xlabel('U(radian)'), ylabel('V(radian)')
%%%%%%%%%% 2

% S=3;    %% scale of downsampling
p=0; 
for m=-200:S:200
    p=p+1; q=0;
    for n=-200:S:200
        q=q+1;
        g(p,q) = f(m+201,n+201);
    end
end

figure 
imshow(g,[])
str_tit = sprintf(' picture after downsampling with S=%d',S);
title(str_tit);

gg= fft2(g);
step=2*pi/length(gg);
[U,V]=meshgrid(0:step:2*pi-step , 0:step:2*pi-step);
figure
mesh(U,V ,abs(gg))
str_tit = sprintf(' Fourier Transfporm of picture after downsampling with S=%d',S);
title(str_tit);
xlabel('U(radian)'), ylabel('V(radian)')


%%%%%%%%%%%  evaluating the effects of leakage and boxcar window

N1 = 10;
N2 = 400;
w = window(@rectwin,N1);
w1 = window(@rectwin,N2); 
wvtool(w1)
wvtool(w)