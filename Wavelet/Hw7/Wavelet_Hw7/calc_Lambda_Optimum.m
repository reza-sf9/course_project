function Lambda_Optimum = calc_Lambda_Optimum(d,level)
L_d = length(d);  %% length d
std_d=std(d);


index=1;
lambda_d=0.001; lambda_u=3*std_d; step=.001;

SURE=zeros(1, length((lambda_u-lambda_d)/step));
for lambda=lambda_d:step:lambda_u
    
    temp1=0; temp2=0;
    
    for i=1:L_d
        
       if abs(d(i))<=  lambda
           temp1=temp1+1;
       end    
       temp2=temp2+(min(abs(d(i)),lambda))^2;
    end
    
    SURE(index) = L_d - 2*temp1 + temp2;

    index=index+1;
end 

t_lambda = lambda_d: step :lambda_u; 
% 
% figure;
% plot(t_lambda,SURE); title(level);

min_SURE=min(SURE);
aa=find(SURE==min_SURE);
Lambda_Optimum= t_lambda(aa(1));  %% the first element of min lambda

end