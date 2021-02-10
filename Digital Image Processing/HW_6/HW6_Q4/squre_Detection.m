function squre_Detection_Hough(f)

f = edge(f);
figure,imshow(f),title('The Edge of Original Image');

[M,N] = size(f);
[ind_x,ind_y] = find(f); %%% find indices that are equal to 1


%%%%%% Calculate Hogh Transform 
step_theta = pi/180;
theta = 0:step_theta:pi-step_theta;
l_theta = length(theta);
M_2 = floor(M*sqrt(2));
ro = zeros(2*M_2 , l_theta);

for i=1:length(ind_x)
    
    temp_ro = round(ind_x(i)* cos(theta) + ind_y(i).*sin(theta));
    ff = zeros(2*floor(M*sqrt(2)) , l_theta);
    for j=1:1:length(temp_ro)
        ff(temp_ro(j)+ M_2,j) = 1;
    end
    ro = ro + ff;
end

%%%%%%% Plot Hough Transform
[NN,MM] = meshgrid(0:179,-M_2:M_2-1);
figure,mesh(MM,NN,ro),xlabel('\rho','FontSize', 20) , ylabel('\theta','FontSize', 20);
title('Hough Transform Of Image')


%%%%%% Extract Maximum Value of Hough and Delete Noisy Peaks
ro_sort = sort(ro(:),'descend');
max_ro = ro_sort(1:10);

p=0;
for i=1:length(max_ro)
    [mm,nn] = find(ro == max_ro(i));
    for k=1:length(mm)
        
        if i==1 && k==1 && abs(nn(k)-180)>5 %%% Ignore theta = 180
            p=p+1;
            ind_x_max(p) =  mm(k);
            ind_y_max(p) = nn(k);
        else
            counter =0;
            for j=1:p
                dd = sqrt((ind_x_max(j)-mm(k)).^2 + (ind_y_max(j)-nn(k)).^2);
                if dd>15
                    counter = counter +1;
                end
            end
            if counter == p && abs(nn(k)-180)>5
                p=p+1;
                ind_x_max(p) =  mm(k);
                ind_y_max(p) = nn(k);
            end
        end
    end
    
end

%%%%% each 2 poits of a square must have the similar Theta
for i=1:length(ind_y_max)
    for j=1:length(ind_y_max)
        if abs(ind_y_max(i) - ind_y_max(j)) < 3
            ind_y_max(i) = ind_y_max(j);
        end
    end
end

%%%%% Delete Noisy Peaks
p=0;
while length(ind_x_max) > 0
    aa = find(ind_y_max(1) == ind_y_max);
    if length(aa) > 1
        for i=length(aa) : -1 : 1
            p=p+1;
            rho_max(p) = round(ind_x_max(aa(i))-M*sqrt(2));
            theta_max(p) = ind_y_max(aa(i))-1;
            ind_x_max(aa(i)) =[];
            ind_y_max(aa(i)) =[];
        end
    else
        ind_x_max(1) = [];
        ind_y_max(1) = [];
    end
end

%%%% Plot Line Between each 2 point with RHO and THETA
gg = zeros(M,N,length(rho_max));
ggg = zeros(M,N);
ggg_1 = zeros(M,N);
figure,
if abs(abs(theta_max(3)-theta_max(1))-90) < 3
    for i=1:length(rho_max)
        
        for m=1:M
            for n=1:N
                if abs(rho_max(i) -( m*cosd(theta_max(i)) + n*sind(theta_max(i)))) < 1
                    gg(m,n,i) = 1;
                    ggg(m,n) = 1;
                    ggg_1(m,n) =1;
                end
            end
        end
        L_gg(i)  = length(find(ggg));
        subplot(2,2,i);
        imshow(ggg),title(['Line with \rho = ',num2str(rho_max(i)),'and \theta = ',num2str(theta_max(i))])
        ggg = zeros(M,N);
    end
end
figure,imshow(ggg_1),title('The result of intersection of 4 Lines');

G  = zeros(M,N,length(rho_max));
for i=1:length(L_gg)
    aa = find(L_gg == max(L_gg));
    G(:,:,i)= gg(:,:,i);
    L_gg(aa) = 0;
end



%%%%%%%%%% Extract 4 Points that obtained from square of 4 Lines
p=0;
for i=1:length(rho_max)
    [m_1 n_1] = find(G(:,:,i)==1);
    for j=i+1: length(rho_max)
        [m_2 n_2] = find(G(:,:,j)==1);
        
        k=0;
        for r=1:length(m_1)
            temp_ind = find(m_1(r) == m_2);
            if length(temp_ind) ~= 0
                for j=1:length(temp_ind)
                    if n_1(r) == n_2(temp_ind(j))
                        p=p+1;
                        pixel(p,:) =[m_1(r) n_1(r)]; %%%%% Cordinate of Corner Points
                        k=1;
                        break
                    end
                end
                if k==1
                    break
                end
            end
        end
        
    end
end

%%%%%% Chack which 2 points can connect to each other
p=0;
for i=1:length(pixel)
    point_1 = pixel(i,:);
    for j=i+1:length(pixel)
        point_2 = pixel(j,:);
        p=p+1;
        mm(p) = round(sqrt( (point_2(2) - point_1(2)).^2  + (point_2(1) - point_1(1)).^2));
        data(p,:) = [point_1 point_2 mm(p)]; 
    end
end

%%%%%%%%%%%%%%% connect 2 points ans make square
g_1=zeros(M,N);
min_d = min(mm);
for i=1:6
    if abs(data(i,5)-min_d)<5
        % distances according to both axes
        y1 = data(i,1);  y2 = data(i,3); x1 = data(i,2);  x2 = data(i,4); 
        xn = abs(x2-x1);
        yn = abs(y2-y1);
        
        % interpolate against axis with greater distance between points;
        % this guarantees statement in the under the first point!
        if (xn > yn)
            xc = x1 : sign(x2-x1) : x2;
            yc = round( interp1([x1 x2], [y1 y2], xc, 'linear') );
        else
            yc = y1 : sign(y2-y1) : y2;
            xc = round( interp1([y1 y2], [x1 x2], yc, 'linear') );
        end
        
        % 2-D indexes of line are saved in (xc, yc), and
        % 1-D indexes are calculated here:
        ind = sub2ind( size(g_1), yc, xc );
        
        % draw line on the image (change value of '255' to one that you need)
        g_1(ind) = 1;
%         figure,imshow(g_1)
    end
end

figure,imshow(g_1)
tit = sprintf(' Lenght of side of square = %d **** Angel of rotation the square = %d \n Cordinate of Corners:\n [m1 n1] = %d %d     [m2 n2] = %d %d   [m3 n3] = %d %d    [m4 n4] = %d %d'...
    ,min_d,min(abs(theta_max)),pixel(1,1),pixel(1,2),pixel(2,1),pixel(2,2),pixel(3,1),pixel(3,2),pixel(4,1),pixel(4,2));
title(tit)

end