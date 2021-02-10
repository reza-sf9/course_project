clc;
clear;
close all;

f = imread('baby2.jpg');
f = rgb2gray(f);

figure(1),imshow(f),title('Original Image');

[M,N] = size(f);
theta = degtorad(38);

A = [cos(theta) -sin(theta) 0; sin(theta) cos(theta) 0; 0  0 1];
%%%%%%%%%%%%%%%%%%%%%%%% offline calculation  size of new matrix(after rotation)
p=0;
for m=1:M
    for n=1:N
        p=p+1;
        U = [m n 1];
        uv(:,p) = (U*A).';
    end
end
u = uv(1,:); v=uv(2,:);
max_u = round(max(u(:))); max_v = round(max(v(:)));
min_u = round(min(u(:))); min_v = round(min(v(:)));

%%%%%%%%%%%%%%%%%%%%%%%%%%%% calculate rotation figure
p=0;
for m=1:M
    for n=1:N
        
        U = [m n 1];
        uv = (U*A).';
        u=round(uv(1,1)); v= round(uv(2,1));
        
        f_rotation(u-min_u+1 , v-min_v+1) = f(m,n);
    end
end

figure(2), imshow(f_rotation),title('Image After Rotation');

%%%%%%%%%%%%%
[x,y] = getpts(figure(2));
X = [x  y ones(6,1)];

[u,v] = getpts(figure(1));
U = [u  v ones(6,1)];

A_Retrieve = inv(X.'*X)*X.'*U;

%%%%%%%%%%%% Retrive Original Image with A Coeficient
%%%%%%%%% get min max
[M1 , N1] = size(f_rotation);
inv_A = inv(A_Retrieve);
uv=[];
p=0;
for m=1:M1
    for n=1:N1
        p=p+1;
        U = [m n 1];
        uv(:,p) = (U*inv_A).';
    end
end
u = uv(1,:); v=uv(2,:);
max_u = round(max(u(:))); max_v = round(max(v(:)));
min_u = round(min(u(:))); min_v = round(min(v(:)));

%%%%%%%%%%%%%%%%%%%%

for m=1:M1
    for n=1:N1
        U = [m n 1];
        uv = (U*inv_A).';
        u=round(uv(1,1)); v= round(uv(2,1));
        
        f_retrieve(u-min_u+2 , v-min_v+1) = f_rotation(m,n);        
    end
end

figure(3), imshow(f_retrieve);
title('retrieve image')

%%%%%%%%%%  crop
[r,c] = find(f_retrieve~=0);
min_r = min(r) + 15; 
min_c = min(c) + 15; 
temp_1 = f_retrieve(min_r,min_c:end);
aa=max(find(temp_1 ~= 0));
max_c  = max(find(temp_1 ~= 0)) + min_c;


temp_2 = f_retrieve(min_r:end,min_c);
max_r = max(find(temp_2 ~= 0)) + min_r;

f_crop = f_retrieve(min_r:max_r , min_c:max_c);
figure(4),imshow(f_crop)
title('cropped retrieve image')


%% Improved Image
%%%%%%%%%%%  Clear Black dots
%%%%%%%%%%%%%%%%%%%%%%%%%% using 8 points around black dots to removo them
[R,C] = find(f_rotation ==0);
%%%% remove indice that row is minimum or moximum on them
remove_min_max_R = find(R~=min(R) & R~=max(R)); 
R = R(remove_min_max_R); C = C(remove_min_max_R);
%%%% remove indice that coloum is minimum or moximum on them
remove_min_max_C = find(C~=min(C) & C~=max(C));
R = R(remove_min_max_C); C = C(remove_min_max_C);

zero_index = [R.';C.'];

f_rotation_Improved = f_rotation;
for i=1:length(zero_index);
    r = zero_index(1,i);
    c = zero_index(2,i);
    temp = f_rotation(r-1:r+1,c-1:c+1);
    index_non_0 = find(temp~=0);
    l_non_0 = length(index_non_0);
    if l_non_0 > 0
           f_rotation_Improved(r,c) = sum(temp(:))/l_non_0;
    end
    
end
 figure(5), imshow(f_rotation_Improved)
 title('improved image after rotation')


%%%%%%
[x,y] = getpts(figure(5));
X = [x  y ones(6,1)];

[u,v] = getpts(figure(1));
U = [u  v ones(6,1)];

A_Retrieve_Improved = inv(X.'*X)*X.'*U;

%%%%%%%%%%%% Retrive Original Image with A Coeficient
%%%%%%%%% get min max
[M1 , N1] = size(f_rotation_Improved);
inv_A = inv(A_Retrieve_Improved);
uv=[];
p=0;
for m=1:M1
    for n=1:N1
        p=p+1;
        U = [m n 1];
        uv(:,p) = (U*inv_A).';
    end
end
u = uv(1,:); v=uv(2,:);
max_u = round(max(u(:))); max_v = round(max(v(:)));
min_u = round(min(u(:))); min_v = round(min(v(:)));

%%%%%%%%%%%%%%%%%%%%

for m=1:M1
    for n=1:N1
        U = [m n 1];
        uv = (U*inv_A).';
        u=round(uv(1,1)); v= round(uv(2,1));
        
        f_retrieve_Improved(u-min_u+2 , v-min_v+1) = f_rotation_Improved(m,n);        
    end
end

figure(6), imshow(f_retrieve_Improved);
title('retrieve improved image')

%%%%%%%%%%%  Clear Black dots
%%%%%%%%%%%%%%%%%%%%%%%%%% using 8 points around black dots to removo them
R=[]; C=[]; remove_min_max_R=[]; remove_min_max_C=[]; zero_index=[];
[R,C] = find(f_retrieve_Improved ==0);
%%%% remove indice that row is minimum or moximum on them
remove_min_max_R = find(R~=min(R) & R~=max(R)); 
R = R(remove_min_max_R); C = C(remove_min_max_R);
%%%% remove indice that coloum is minimum or moximum on them
remove_min_max_C = find(C~=min(C) & C~=max(C));
R = R(remove_min_max_C); C = C(remove_min_max_C);

zero_index = [R.';C.'];

f_retrieve_Improved2 = f_retrieve_Improved;
for i=1:length(zero_index);
    r = zero_index(1,i);
    c = zero_index(2,i);
    temp = f_retrieve_Improved(r-1:r+1,c-1:c+1);
    index_non_0 = find(temp~=0);
    l_non_0 = length(index_non_0);
    if l_non_0 > 0
           f_retrieve_Improved2(r,c) = sum(temp(:))/l_non_0;
    end
    
end
 figure(7), imshow(f_retrieve_Improved2)
 title('removing black dots from retrieved image')



%%%%%%%%%%%  crop
[r,c] = find(f_retrieve_Improved2~=0);
min_r = min(r)+15; 
min_c = min(c)+15; 
temp_1 = f_retrieve_Improved2(min_r,min_c:end);
max_c  = max(find(temp_1 ~= 0)) + min_c;


temp_2 = f_retrieve_Improved2(min_r:end,min_c);
max_r = max(find(temp_2 ~= 0)) + min_r;

f_crop = f_retrieve_Improved2(min_r:max_r , min_c:max_c);
figure(8),imshow(f_crop)
title('croped retrieve image')

error_1 = norm(A(1:2,1:2)-A_Retrieve(1:2,1:2))*100; %%% norm of error without removing black dots 
error_2 = norm(A(1:2,1:2)-A_Retrieve_Improved(1:2,1:2))*100; %%% norm of error after removing black dots
