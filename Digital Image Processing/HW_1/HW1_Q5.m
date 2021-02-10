clc;
clear ;
close all



A1 = ones(10,5);
A2 = 2*ones(10,5);

X = [A1 A2];                  %% input signal
H1 = [1 0 -1;1 0 -1;1 0 -1];   %% filrers
H2 = [1 1 1 ; 0 0 0 ; -1 -1 -1];

X=X;  %%%% chose input
H=H1;   %%%% chose filter 
l_m = length(X) + length(H);
zeros(l_m,l_m);

figure, imshow(H1,[]),title('h1')
figure, imshow(H2,[]),title('h2');

m=-1;
n=-1;
y(m+2,n+2)=0;
% 
for m= -floor(length(H)/2):length(X)-1 +floor(length(H)/2)
    for n= -floor(length(H)/2):length(X)-1 +floor(length(H)/2)
        row_index = n + ceil(length(H)/2);
        col_index = m + ceil(length(H)/2);
        y(row_index, col_index) =0;
        for k1=0:length(X)-1
            for k2=0:length(X)-1
                a=n-k1 + ceil(length(H)/2);
                b=m-k2 + ceil(length(H)/2);
                if    min( a , b) >= 1   &&   max(a,b)  <= length(H)
                    y(row_index, col_index) =  y(row_index, col_index) + X(k1+1,k2+1)* H(a, b);
                end
            end
        end
    end
end
y
figure, imshow(X,[])
figure,imshow(y,[]), title('h1**x''')

