function Kass_Method(f , y , x, alpha , beta , gamma , win_length)
a=alpha;
b=beta;
c=gamma;
w_2 = floor(win_length/2);

%%%% filetr Image with a Gaussian
l_win = 15; sigma = l_win/6;
h=fspecial('gaussian', l_win , sigma);
f=imfilter(f,h);


N=length(x);
aa = beta;
%%%%% Coeficients to caculate E
alpha=alpha*ones(1,N);
beta=beta*ones(1,N);
beta2=aa*ones(1,N);
gamma=gamma*ones(1,N);

%%%%%%% Calculate Gradient of Image
[grad_x , grad_y] =  gradient(f);
Grad = grad_x.^2 + grad_y.^2;

%%%%% Threshold for corner detection citeria
max_G = max(Grad(:)); min_G = min(Grad(:)); d = max_G - min_G;
grad_thr = (.7*max(Grad(:)) - .4*min(Grad(:)) - min_G)/d;

curv_thr = .5;

%%%%%%%%%%%%%%%%
x_1 = x;
y_1 = y;
count_change = 1;
while  count_change > 0
    
    %     while(1)
    %%%%%%%%%%%%%%%%%%%%%%%%
    
    x_aug = [x_1(N) x_1.' x_1(1)];
    y_aug = [y_1(N) y_1.' y_1(1)];
    d = 0;
    p1=0;
    for i=2:N+1
        p1=p1+1;
        d = d + sqrt((x_aug(i) - x_aug(i-1)).^2 + (y_aug(i) - y_aug(i-1)).^2);
    end
    d = d/N;
    
    count_change = 0;
    for i=2:N+1
        x_aug = [x_1(N) x_1.' x_1(1)];
        y_aug = [y_1(N) y_1.' y_1(1)];
        
        
        for m=-w_2:w_2
            for n=-w_2:w_2
                X = x_aug(i) + m;
                Y = y_aug(i) + n;
                E_Count(m+w_2+1,n+w_2+1) = (d-sqrt((X-x_aug(i-1)).^2 + (Y-y_aug(i-1)).^2)).^2;
                E_Curv(m+w_2+1,n+w_2+1) = (x_aug(i+1) - 2*X + x_aug(i-1)).^2 + (y_aug(i+1) - 2*Y + y_aug(i-1)).^2;
            end
        end
        E_Grad= Grad(x_aug(i)-w_2:x_aug(i)+w_2 , y_aug(i)-w_2:y_aug(i)+w_2);
        
        %%%% NOrmlized
        E_Count_norm = normlize(E_Count,1);
        E_Curv_norm  = normlize(E_Curv,1);
        E_Grad_norm = normlize(E_Grad,2);
        
        
        
        temp = alpha(1,i-1).*E_Count_norm + beta(1,i-1).*E_Curv_norm - gamma(1,i-1).*E_Grad_norm;
        
        
        
        [xx,yy] = (find(temp == min(temp(:))));
        x_n = x_1(i-1) + xx - (w_2+1);
        y_n = y_1(i-1) + yy - (w_2+1);
        
        if x_n ~= x_1(i-1) || y_n~=y_1(i-1)
            count_change = count_change + 1;
        end
        
        %%%%%%%%%% Update Points
        x_1(i-1) = x_n;
        y_1(i-1) = y_n;
                E_Curv_final(i-1) = E_Curv_norm(xx , yy); 
                E_Grad_final(i-1) = E_Curv_norm(xx , yy);
    end
    
    E_Curv_final = [ E_Curv_final(end) E_Curv_final E_Curv_final(1) ];
    
    for i=2:N+1
        if E_Curv_final(i) > E_Curv_final(i+1)  && ...
                E_Curv_final(i) > E_Curv_final(i-1)&&...
                E_Curv_final(i) > curv_thr   && E_Grad_final(i-1)>grad_thr
            
            beta(i)=0;
        end
    end
    
    %%%%%% Plot
    hold off
    imshow(f,[]),hold on,plot([y_1;y_1(1)],[x_1;x_1(1)],'r','LineWidth',3)
    title(['\alpha = ',num2str(a),'  \beta = ',num2str(b),'  \gamma = ',num2str(c),'  win Length = ', num2str(win_length)])
    pause(.1)
end

end