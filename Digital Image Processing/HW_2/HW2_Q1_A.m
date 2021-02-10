clc;
clear;
close all;

%%%%%%%%%%%%%%%%%%% Matrix A Interpolation
A =[-0.5 1.5 -1.5 0.5; 1 -2.5 2 -0.5;-0.5 0 0.5 0;0 1 0 0];
f = [0 0 1 3 -2 5 5 5 -1 3 3 2 0 0];

step=1/4.5;

p_n = (3:-1:0).';
p=0;
for n=2:step:length(f)-2
    p = p+1;
    f_n(p) = fix(n);               %%% fixed value of n
    d_n = n-f_n(p);                %%% delta_n
    points = f(f_n(p)-1:f_n(p)+2); %%% extract 4  points in order to interolate
    Coef = A*points.';             %%% calculate coeficients
    f_A(p) = Coef.' * d_n.^ p_n;   %%% calulate interpolated point
    
end

ind_5=find(f_n==5);                          %%% indices include f(5)
ind_7=find(f_n==7);                          %%% indices include f(7)
f_A_5_7 = f_A(ind_5(1):ind_7(end))         %%% interpolated Data Between 5 and 7



%%%%%%%%%%%%%%%% Cubic Spline Interpolation 

Cubic_Coef = load('Cubic_Coef.txt')         %% Coeficient for Cubic Spline Method

%%% Y is Si(xi)
for j=1:12
    Y(j) = f(j) - 2*f(j+1) + f(j+2);
end

M = 6.* inv(Cubic_Coef)*Y.' ;                %% M is the second derivative matrice [Si"(xi)]
M = [M(1);M;M(12)];                          %% Update M and calculate M_1 and M_14

%%%%% calculate a , b , c , d 
M_1 = M(2:end);
f_1 = f(2:end);

M_0 = M(1:12); M_1 = M_1(1:12); 
f_0 = f(1:12).'; f_1 = f_1(1:12).';

a = (M_1 - M_0)/6;
b = M_0/2;
c = (f_1 - f_0) - ((M_1 + 2*M_0)/6);
d = f_0;


p=0;
for n=2:step:length(f)-2
    p = p+1;
    n_Cubic(p) = fix(n);
    d_n = n-n_Cubic(p);
    points = f(n_Cubic(p)-1:n_Cubic(p)+2);
    Coef = [a(n_Cubic(p)) ;b(n_Cubic(p)) ;c(n_Cubic(p)) ;d(n_Cubic(p))];   
    f_Cubic(p) = Coef.' * d_n.^ p_n;
end

ind_5=find(n_Cubic==5);                     %%% indices include f(5)
ind_7=find(n_Cubic==7);
%%% indices include f(7)
f_Cubic_5_7 = f_Cubic(ind_5(1):ind_7(end)) %%% interpolated Data Between 5 and 7

%%%% Calculate rms error 
rms_Error = norm(f_Cubic - f_A)/norm(f_Cubic)*100

%%%% Plotiing Cubic and Aproximated Method
t = (2:step:length(f)-2);
figure;
plot(t,f_Cubic,'r'); hold on;
plot(t,f_A,'b--'); hold on;
plot(2:length(f)-2,f(2:length(f)-2),'k.','MarkerSize',30);
xlabel('n'), ylabel('Interplated Line') , title('Comparision Cubic Spline & Aproximated Spline Method')
legend('Cubic Spline','A matrix Spline','Sequence')


