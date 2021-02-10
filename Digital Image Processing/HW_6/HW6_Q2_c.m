clc;
clear;
close all;


f = imread('bloodcel_95.jpg');
f = rgb2gray(f);
f = im2double(f);
[M,N] = size(f);
f = f + .01.*randn(M,N);

figure,imshow(f), title('Original Image')

%%%%%% SUSSAN
for w = 11:1:11
    for thr_sussan =.105:.005:.105
        
        N_Length = 2*w + 1;
        
        %%%% zero-pad f
        ff = zeros(M+2*w , N+2*w);
        ff(w+1:end-w , w+1:end-w) = f;
        
        
        Nucleus = zeros(N_Length , N_Length);
        
        for m=-w:w
            for n=-w:w
                r = sqrt(m.^2 + n.^2);
                if ceil(r) <= w
                    Nucleus(m+w+1,n+w+1) = 1;
                end
            end
        end
         figure, imshow(Nucleus);
        
        n_max = length(find(Nucleus));
        p =0;
        for m=w+1:M-(w+1)
            p=p+1; q=0;
            for n=w+1:N-(w+1)
                q=q+1;
                center_pixel = ff(m,n);
                
                temp = zeros(N_Length,N_Length);
                for mm = -w:w
                    for nn = -w:w
                        temp(mm+w+1,nn+w+1) = Nucleus(mm+w+1,nn+w+1) .*  abs(center_pixel - ff(m+mm,n+nn));
                    end
                end
                temp_2(p,q) = length(find(temp < thr_sussan));
                
                
            end
        end
        temp_3 = .75*n_max - temp_2;
        ind = find(temp_2 < .75*n_max);
        [MM,NN] = size(temp_2);
        temp_4 = zeros(MM,NN);
        temp_4(ind) = temp_3(ind);
        
        fig=figure;
        imshow(temp_3),title(['Neucles , length =',num2str(N_Length),'  thr = ',num2str(thr_sussan)]);
        name = sprintf('sussan_Length_%d__Thr_%s.jpg',N_Length,num2str(thr_sussan));
        saveas(fig,name);
        g = temp_4;
        save('f_WLennght_23_Thr_.105.mat','g');
    end
end
