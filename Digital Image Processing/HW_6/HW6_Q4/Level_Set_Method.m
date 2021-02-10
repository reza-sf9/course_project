function Level_Set_Method(f)

[M N]=size(f);
[fx,fy]=gradient(f);
absgrad=sqrt(fx.^2+fy.^2);
gama=2;
g_m=1+gama*absgrad;
p=1./g_m;
p=double(p);
phi=ones(M,N);
for m=1:M
    for n=1:N
        D=sqrt((m-100)^2+(n-100)^2);
        if D<=75,phi(m,n)=-1;
        end
    end
end

phi = double((phi > 0).*(bwdist(phi < 0)-0.5) - (phi < 0).*(bwdist(phi > 0)-0.5));
figure;surf(phi);

e=.51;
figure,
for t=1:77
phi=phi+e.*p;    
    phi = double((phi > 0).*(bwdist(phi < 0)-0.5) - (phi < 0).*(bwdist(phi > 0)-0.5));
    
    imshow(f,[]);hold on
    z=contour(phi,[0,0],'r','LineWidth',3);title(['number of iterations=',num2str(t)]);
    pause(0.00001)
end


end