%% In the name of God

%% DIP - HW03 - 9400365 - Mohammad Hossein Bahonar

clear
clc
close all

f = imread('text_10.jpg');                                  % Reading the image
f = rgb2gray(f);

%%%% Matlab OTSU Method
level = graythresh(f);
mat_otsu = level*256;
BW = im2bw(f,level);
figure,imshow(BW),title('Matlab OTSU Method')


f = 255 * im2double(f);



WinSize = 25;                                               % Window Size
K_Sauvola = 0.5;                                            % Sauvola coefficient
K_Niblack = -0.2;                                           % Niblack coefficient

[M,N] = size(f);
fPreIntImage  = zeros(M+1,N+1);                             % Zero padded image
fPreIntImage(2:M+1,2:N+1) = f;
f2PreIntImage = fPreIntImage.^2;                            % Zero padded squared image

fIntImage  = zeros(M+1,N+1);                                % Initialization of integral sum of F
f2IntImage = fIntImage;                                     % Initialization of integral sum of F 
for m=2:M+1                                                 %  squared
    for n=2:N+1
        fIntImage(m,n)  = fIntImage (m,n-1)   + fIntImage    (m-1,n)...     % Calculating integral sum
                        - fIntImage (m-1,n-1) + fPreIntImage (m,n);         %  of the image
        f2IntImage(m,n) = f2IntImage(m,n-1)   + f2IntImage   (m-1,n)...     % Calculating integral sum
                        - f2IntImage(m-1,n-1) + f2PreIntImage(m,n);         %  of image squared
    end
end
fIntImage  = fIntImage (2:M+1,2:N+1);
f2IntImage = f2IntImage(2:M+1,2:N+1);

L = fix(WinSize/2);

%% Sauvola Section
fSauvola = zeros(M,N);                                  % Sauvolva initialization

for m=L+1:M-L
     for n=L+1:N-L
         
         Mean = fIntImage(m+L,n+L)  - fIntImage(m+L,n-L)...     % Calculating mean using integral sum
               -fIntImage(m-L,n+L)  + fIntImage(m-L,n-L);
         Mean = Mean / (WinSize^2);                             % Normalizing Mean
         Std  = f2IntImage(m+L,n+L) - f2IntImage(m+L,n-L)...    % Calculating std using integral sum
               -f2IntImage(m-L,n+L) + f2IntImage(m-L,n-L);
         Std_zig = Std / (WinSize-1)^2 - Mean^2;                % Normalizing std
         Std_zig = sqrt(Std_zig);
         Threshold = Mean * (1 - K_Sauvola*(1-Std_zig/127) );   % Calculating Sauvola threshold
         if f(m,n)<=Threshold                                   % Thresholding image
             fSauvola(m,n)=0;
         else
             fSauvola(m,n)=1;
         end
         
     end
end
 
% figure, imshow(fSauvola,[])                                     % Displaying image
% title(['Sauvola thresholding, WinSize=' num2str(WinSize) ', K=' num2str(K_Sauvola)])

%% Niblack Section
fNiblack = zeros(M,N);                                          % Niblack initialization

for m=L+1:M-L
     for n=L+1:N-L
         
         Mean = fIntImage(m+L,n+L)  - fIntImage(m+L,n-L)...     % Calculating mean using integral sum
               -fIntImage(m-L,n+L)  + fIntImage(m-L,n-L);
         Mean = Mean / (WinSize^2);                             % Normalizing mean
         Std  = f2IntImage(m+L,n+L) - f2IntImage(m+L,n-L)...    % Calculating std using integral sum
               -f2IntImage(m-L,n+L) + f2IntImage(m-L,n-L);
         Std_zig = Std / (WinSize-1)^2 - Mean^2;                % Normalizing std
         Std_zig = sqrt(Std_zig);
         Threshold = Mean + K_Niblack * Std_zig;                % Calculating threshold of Niblack
         if f(m,n)<=Threshold                                   % Thresholding image
             fNiblack(m,n)=0;
         else
             fNiblack(m,n)=1;
         end
         
     end
end
%  
% figure, imshow(fNiblack,[])
% title(['Niblack thresholding, WinSize=' num2str(WinSize) ', K=' num2str(K_Niblack)])

%% Iterative Global Thresholding
[hDist,~] = hist( f(:),0:255 );                                 % Histogram of image

Threshold = 127;                                                % Initial threshold
Finished = 0;
while Finished==0
    
    Mean1 = (0:Threshold-1) * hDist(1:Threshold).';             % Mean of the pixels below threshold
    Mean1 = Mean1 ./ sum( hDist(1:Threshold) );
    Mean2 = (Threshold:255) * hDist(Threshold+1:256).';         % Mean of the pixels above threshold
    Mean2 = Mean2 ./ sum( hDist(Threshold+1:256) );
    ThresholdNew = (Mean1+Mean2)/2;                             % Calculating new threshold
    if abs(Threshold-ThresholdNew)<1, Finished=1; end           % Stopping criteria
    Threshold = round(ThresholdNew);
    
end

fGlobalThreshold = f > Threshold;                               % Thresholding

fGlobalThreshold(end,:)=0; fGlobalThreshold(1,:)=0;             % Image displaying
fGlobalThreshold(:,end)=0; fGlobalThreshold(:,1)=0;
% figure, imshow(fGlobalThreshold)
% title(['Thresholding with Iterative Method, T=' num2str(Threshold)])

%% OTSU

[hDist,~] = hist( f(:),0:255 );
SumTotal  = sum(hDist);                                      
Error     = zeros(1,255);                                   % Error related to different thresholds
for ThCurrent=1:255
    SumPart1 = sum( hDist(1:ThCurrent) );
    SumPart2 = sum( hDist(ThCurrent+1:256));
    Prob1    = hDist(1:ThCurrent)     / SumPart1;           % Probability related to first part pixels
    Prob2    = hDist(ThCurrent+1:256) / SumPart2;           % Probability related to second part pixel
    Mean1    = (1:ThCurrent)*Prob1';                        % Mean of pixels below threshold
    Mean2    = (ThCurrent+1:256)*Prob2';                    % Mean of pixels above threshold
    
    %%%%%%% I think this part is wrong
    Std1     = sum( (hDist(1:ThCurrent)-Mean1).^2 ) / SumPart1;             % Std of pixel below thre
    Std2     = sum( (hDist(ThCurrent+1:256)-Mean2).^2) / SumPart2;          % Std of pixel above thre
    Error(ThCurrent) = SumPart1/SumTotal*Std1 +SumPart2/SumTotal*Std2;      % OTSU objective function
    
end
figure, plot(Error)
title('Error Of OTSU Method')
xlabel('Pixel Value')
ylabel('Error Probability')
[~,IndMin] = min(Error);                                    % Optimum OTSU threshold

fOTSU = f>IndMin;                                           % OTSU Thresholding
fOTSU(1,:)=0;
fOTSU(end,:)=0;
fOTSU(:,1)=0;
fOTSU(:,end)=0;
figure, imshow(fOTSU,[])
title(['Output of OTSU Method, T=' num2str(IndMin)])


k
