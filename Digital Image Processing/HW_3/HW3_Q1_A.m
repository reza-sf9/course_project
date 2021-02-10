clc;
clearvars -except fIntImage
close all;

I=imread('degarded_text.png');
I = 255 * im2double(I);                %%% conver unit8 to decimal number

figure, imshow(I,[]),title('Original Image');

%%%%%%%%%%%%%%%%%%%%%%%%%%     NiBlack & Savoula Approach   %%%%%%%%%%%%%%
k=.2;
R = 74;
win_length = 7;

[M,N] = size(I);
step_win = fix(win_length/2);

I_temp = I; %%%%%%%%%%%%%%%%%%%
[M_temp,N_temp] = size(I_temp);

step_win = fix(win_length/2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Integeral Image  %%%%%%%%%%%%%%%%%%%%%%%%%

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
p=0;
for m=1 + step_win : M - step_win
    p=p+1; q=0;
    for n=1 + step_win :N - step_win
        q=q+1;
        m_1 = m+ step_win+1;  n_1 = n + step_win+1;            coef_1 = m_1*n_1;
        m_2 = m_1 - win_length;   n_2 = n_1;                   coef_2 = m_2*n_2;
        m_3 = m_1;                n_3 = n_1 - win_length;      coef_3 = m_3*n_3;
        m_4 = m_1-win_length;     n_4 = n_1-win_length;        coef_4 = m_4*n_4;
        w_size = coef_1+ coef_4 - coef_2 - coef_3;  %%% window size
        
        mu_Int(p,q) = 1/win_length^2 * (Int_Mat(m_1,n_1) - Int_Mat(m_2,n_2) - Int_Mat(m_1,n_3) +  Int_Mat(m_4,n_4));
        
        var_Int(p,q) = sqrt((1/win_length^2 * (Int_Mat_2(m_1,n_1) - Int_Mat_2(m_2,n_2) - Int_Mat_2(m_1,n_3) +  Int_Mat_2(m_4,n_4)))...
            - mu_Int(p,q).^2);
        
        T_Niblack_Int(p,q) = mu_Int(p,q) - k*var_Int(p,q);            %%% Threshold Niblack
        T_Savoula_Int(p,q) = mu_Int(p,q)*(1+k*(var_Int(p,q)/R - 1));  %%% THreshold Savoula
    end
end

T_temp = zeros(M,N);
T_temp(1 + step_win : M - step_win,1 + step_win :N - step_win) = T_Niblack_Int;
T_Niblack_Int = T_temp;

T_temp = zeros(M,N);
T_temp(1 + step_win : M - step_win,1 + step_win :N - step_win) = T_Savoula_Int;
T_Savoula_Int = T_temp;


%%%% Niblack
I_Niblack=I;
index_Ni_255 = find(I_Niblack >= T_Niblack_Int);
index_Ni_0 =   find(I_Niblack < T_Niblack_Int);

I_Niblack(index_Ni_255)=255;
I_Niblack(index_Ni_0)=0;
temp=I_Niblack(4:end-4,4:end-4);
figure, imshow(I_Niblack);
str_tit =sprintf('Niblack Thresholding , Window Length =%d ,Parametere K =%.1f',win_length,k); title(str_tit);

%%%%%% Savoula
I_Savoula = I;
index_Sv_255 = find(I_Savoula >= T_Savoula_Int);
index_Sv_0 =   find(I_Savoula < T_Savoula_Int);

I_Savoula(index_Sv_255)=255;
I_Savoula(index_Sv_0)=0;
figure, imshow(I_Savoula);
str_tit = sprintf('Savoula Thresholding , Window Length =%d ,Parametere K =%.1f , R =%d',win_length,k,R); title(str_tit);