clc;
clear;
close all;

step =1;
p=0;
for m=0:step:149
    p=p+1 ; q=0;
    for n=0:step:149
        q=q+1;
        
        %         f(p,q) = cos(2*pi*m/7.5)*cos(2*pi*n/10);
        f(p,q) = cos((2*pi/7.5*m));
    end
end
figure;
imshow(f,[])


figure
mesh(f)
xlabel('n'),ylabel('m')

F = fft2(f);
l_f=length(F);
step_f = 2*pi/l_f;
[M_F,N_F]=meshgrid( 0:step_f:2*pi-step_f,  0:step_f:2*pi-step_f);
figure;
mesh(M_F , N_F ,abs(F));
xlabel('n freq (V)'); ylabel('m freq (U)')


%%%%%%%%%%%%%%%%%
theta = -pi/3; T = 7.5;
p=0;
for m=0:step:149
    p=p+1 ; q=0;
    for n=0:step:149
        q=q+1;
        
                 mm=m*cos(theta)-n*sin(theta);
                 f(p,q) = cos((2*pi/T).*mm);
    end
end
figure;
imshow(f,[])


figure
mesh(f)
xlabel('n'),ylabel('m')

[M , N] = size(f);
ff =  zeros(M,N);
for m=1:M
    for n=1:N
        ff(m,n) = ((-1)^(m+n)) * f(m,n);
    end
end




F = fft2(ff);
l_f=length(F);
step_f = 2*pi/l_f;
[M_F,N_F]=meshgrid( 0:step_f:2*pi-step_f,  0:step_f:2*pi-step_f);
figure;
mesh(M_F , N_F ,abs(real(F)));
xlabel('n freq (V)'); ylabel('m freq (U)')
pa = real(abs(F));
[a , b] = find(abs(real(F)) > 4500);


aa= mean((a(1)-75));
bb= mean((b(1)-75));
atand(aa/bb)

j