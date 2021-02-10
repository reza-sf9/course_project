clc
close all
clear

f=imread('degarded_text.png');
f = 255 * im2double(f);
[M,N] = size(f);

r_max = sqrt(220^2 + 150^2);  r_min = 0;
T_max = 100;                  T_min = 10;

p=0; 
for m=-150:149
    p=p+1; q=0;
    for n=-220:219
        q=q+1;
        r= sqrt(m^2 + n^2);
        T = ((T_min -T_max)/(r_max - r_min))*r + T_max;
        chirp(p,q) = sin(2*pi*r/T);
    end
end
chirp=fix(128*(chirp+1));  %%% Normlized noise between 0 - 255

g= .5*(f+chirp);
figure,imshow(g,[]),title('f + noise (chirp)')


%%%%%%%%%%%%%%%%%%%%%%%%%%%     NiBlack & Savoula Approach   %%%%%%%%%%%%%%
I = g;

k=.45;             %%% k for Niblack adn Savoula Method
R = 128;           %%% R for Savoula Method
win_length = 5;   %%% window length
    
[M,N] = size(I);

I_temp = I; 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Integeral Image  %%%%%%%%%%%%%%%%%%%%%%%%%
% for R=20 : 10 : 250
step_win = fix(win_length/2);

I_Integral = zeros(M+1,N+1);
I_Integral(2:M+1,2:N+1) = I;

I_Integral_2 = zeros(M+1,N+1);
I_Integral_2(2:M+1,2:N+1) = I.^2;

Int_Mat = zeros(M+1,N+1);          %%% Integral matrix
Int_Mat_2 = zeros(M+1,N+1);        %%% Integral matrix power 2

for m=2:M+1
    for n=2:N+1
        Int_Mat(m,n) = Int_Mat(m-1,n) + Int_Mat(m,n-1) - Int_Mat(m-1,n-1) +  I_Integral(m,n);
        Int_Mat_2(m,n) = Int_Mat_2(m-1,n) + Int_Mat_2(m,n-1) - Int_Mat_2(m-1,n-1) +  I_Integral_2(m,n);
    end
end

%%% calculate mean and variance
p=0; T_Niblack =[]; T_Savoula = [];
for m=1 + step_win : M - step_win
    p=p+1; q=0;
    for n=1 + step_win :N - step_win
        q=q+1;
        m_1 = m+ step_win+1;  n_1 = n + step_win+1;    coef_1 = m_1*n_1;
        m_2 = m_1 - win_length;   n_2 = n_1;                   coef_2 = m_2*n_2;
        m_3 = m_1;                n_3 = n_1 - win_length;      coef_3 = m_3*n_3;
        m_4 = m_1-win_length;     n_4 = n_1-win_length;        coef_4 = m_4*n_4;
        mu = 1/win_length^2 * (Int_Mat(m_1,n_1) - Int_Mat(m_2,n_2) - Int_Mat(m_1,n_3) +  Int_Mat(m_4,n_4));
        
        var = sqrt((1/win_length^2 * (Int_Mat_2(m_1,n_1) - Int_Mat_2(m_2,n_2) - Int_Mat_2(m_1,n_3) +  Int_Mat_2(m_4,n_4)))...
            - mu.^2);
        
        T_Niblack(p,q) = mu - k*var;            %%% Threshold Niblack
        T_Savoula(p,q) = mu*(1+k*(var/R - 1));  %%% THreshold Savoula
    end
end

T_temp = zeros(M,N);
T_temp(1 + step_win : M - step_win ,1 + step_win :N - step_win ) = T_Niblack;
T_Niblack = T_temp;

T_temp = zeros(M,N);
T_temp(1 + step_win : M - step_win,1 + step_win :N - step_win) = T_Savoula;
T_Savoula = T_temp;


%% Niblack
I_Niblack=I;
index_Ni_255 = find(I_Niblack >= T_Niblack);
index_Ni_0 =   find(I_Niblack < T_Niblack);

I_Niblack(index_Ni_255)=255;
I_Niblack(index_Ni_0)=0;
temp=I_Niblack(4:end-4,4:end-4);
figure, imshow(I_Niblack);
str_tit =sprintf('Niblack Thresholding , Window Length =%d ,Parametere K =%.2f',win_length,k); title(str_tit);

%%%%%%% Savoula
I_Savoula = I;
index_Sv_255 = find(I_Savoula >= T_Savoula);
index_Sv_0 =   find(I_Savoula < T_Savoula);

I_Savoula(index_Sv_255)=255;
I_Savoula(index_Sv_0)=0;
figure, imshow(I_Savoula);
str_tit = sprintf('Savoula Thresholding , Window Length =%d ,Parametere K =%.2f , R=%d',win_length,k,R); title(str_tit);


%%%% save as 
% str_name = sprintf('R_%d',R);
% str_name = sprintf('Savoula_R_%.2f.jpg',R);
% saveas(h,str_name,'jpg')
% end