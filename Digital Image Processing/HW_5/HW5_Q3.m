clc
clear
close all;

load('p5_1.mat');
figure, imshow(f),title('Original Image')

[M,N]=size(f);

ff =  zeros(M,N);
for m=1:M
    for n=1:N
        ff(m,n) = ((-1)^(m+n)) * f(m,n);
    end
end


F = fft2(ff);
step_u = 2*pi/M; step_v = 2*pi/N;
[M_F,N_F]=meshgrid( -pi:step_v:pi-step_v,  -pi:step_u:pi-step_u);
figure;
mesh(M_F , N_F ,abs(real(F)));
xlabel('n freq (V)'); ylabel('m freq (U)'),title('Fourier Transform of Image')


%%%%%% Find max 
G= F;
f_m = floor(M/2); f_n = floor(N/2);
G(f_m-20 : f_m+20 , f_n -20:f_n+20) =0; %%%% remove the frequency content areound zero 
GG = sort(abs(real(G(:))),'descend');
G_max_ind= length(find(GG> .5*max(GG)));
GG = GG(1:G_max_ind);
for i=1: 2 :length(GG)
    [temp_1 , temp_2] = find(abs(real(F)) == GG(i));
    ind_m(i)= temp_1(1);      ind_n(i)= temp_2(1); 
    ind_m(i+1)= temp_1(2);    ind_n(i+1)= temp_2(2); 
end

ind_1_4 = find(ind_n > N/2);
ind_m_2 = ind_m(ind_1_4) - M/2;
ind_n_2 = ind_n(ind_1_4) - N/2;

for i=1 : length(ind_1_4)
    thet(i) = atand(ind_m_2(i) / ind_n_2(i));
    omega_m(i) = abs(ind_m_2(i)*pi/(M/2));
    omega_n(i) = abs(ind_n_2(i)*pi/(N/2));
end


%% %%%%%%%%%%%    Ideal Band Reject Filter
D0=[]; p=0; q=0;
W= 10; %%% 10 pixel
for i=1:length(ind_m)
    dist = sqrt((ind_m(i)-M/2)^2 + (ind_n(i)-N/2)^2);
    count = 0;
    %%% check this distanse existed before or not
    [r,c] = size(D0);
    for j=1: r
        dif = D0(j,1)-dist;
        if abs(dif) <= W
            q=q+1;
            D0(j,1) = (D0(j,1) + dist)/2;
            count = count +1;
            ind_pa(1,1,q) = ind_m(D0(j,2));  ind_pa(1,2,q) = ind_n(D0(j,2));
            ind_pa(2,1,q) = ind_m(i);        ind_pa(2,2,q) = ind_n(i);
        end
    end
    
    if count ==0
        p=p+1;
        D0(p,1) = dist;
        D0(p,2) = i;     %%%index number
    end
    
end
D0 = D0(:,1);


%%%%%%%%%%% Ideal Band Reject
for i=1:length(D0)
    H = zeros(M , N);
    for u=1:M
        for v=1:N
            D = sqrt((u-M/2)^2 + (v-N/2)^2);
            if D <= D0(i)-W/2 || D >= D0(i)+W/2
                H(u,v) = 1;
            end
        end
    end
    H_Ideal(:,:,i) = H;
end

 figure,
for i=1:length(D0)
    temp = H_Ideal(:,:,i);
    step_u = 2*pi/M; step_v = 2*pi/N;
    [M_F,N_F]=meshgrid( 0:step_v:2*pi-step_v,  0:step_u:2*pi-step_u);
    subplot(1,3,i);
    mesh(M_F , N_F ,temp);
    xlabel('n freq (V)'); ylabel('m freq (U)'),title(['Ideal Band Reject - filter #', num2str(i)])
end

%%%%%
F_hat= F;
for i=1:length(D0)
    temp = H_Ideal(:,:,i);
    F_hat = F_hat.* temp;
    step_u = 2*pi/M; step_v = 2*pi/N;
    [M_F,N_F]=meshgrid( 0:step_v:2*pi-step_v,  0:step_u:2*pi-step_u);
    figure;
    mesh(M_F , N_F ,abs(real(F_hat)));
    xlabel('n freq (V)'); ylabel('m freq (U)'),title(['Fourier Transform Reesponse after Mutiply to filter #', num2str(i),' ** Ideal Band Reject'])
end



g = real(ifft2(F_hat));
for m=1:M
    for n=1:N
        gg(m,n) = ((-1)^(m+n)) * g(m,n);
    end
end
figure, imshow(gg),title(['denoising Image with Ideal Band Reject - with W = ',num2str(W)])


%%%%%%%%%%% Butterworth Notch filter 
d0=25; n=4;
[r,c,d] = size(ind_pa);
for i=1 : d
    temp_1 = ind_pa(1,:,i);
    temp_2 = ind_pa(2,:,i);
    for u=1:M
        for v=1:N
    D_1 = sqrt( (u-temp_1(1))^2 + (v-temp_1(2))^2);
    D_2 = sqrt( (u-temp_2(1))^2 + (v-temp_2(2))^2);
    H(u,v) = inv(1+ (d0^2/(D_1*D_2))^n);
        end
    end
    H_Butter(:,:,i) = H; 
end


for i=1:d
    temp = H_Butter(:,:,i);
    step_u = 2*pi/M; step_v = 2*pi/N;
    [M_F,N_F]=meshgrid( 0:step_v:2*pi-step_v,  0:step_u:2*pi-step_u);
%     figure,    
%     mesh(M_F , N_F ,temp);
%     xlabel('n freq (V)'); ylabel('m freq (U)'),title(['ButterWorth Notch filter - filter #', num2str(i),' ** D0=',num2str(d0)])
end

%%%%%
F_hat= F;
for i=1:d
    temp = H_Butter(:,:,i);
    F_hat = F_hat.* temp;
    step_u = 2*pi/M; step_v = 2*pi/N;
    [M_F,N_F]=meshgrid( 0:step_v:2*pi-step_v,  0:step_u:2*pi-step_u);
    figure;
    mesh(M_F , N_F ,abs(real(F_hat)));
    xlabel('n freq (V)'); ylabel('m freq (U)'),title(['Fourier Transform Reesponse after Mutiply to filter #', num2str(i),' ** ButterWorth Notch filter'])
end


g = real(ifft2(F_hat));
for m=1:M
    for n=1:N
        gg(m,n) = ((-1)^(m+n)) * g(m,n);
    end
end
figure, imshow(gg),title(['denoising Image with ButterWorth Notch filter - D0 = ',num2str(d0)])










