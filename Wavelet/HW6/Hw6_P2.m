clc;
clear;
close all;

%%%%%%%%%%%%%%%%%%%%% Problem 2  %%%%%%%%%%%%%%%%

x=zeros(1,20);
index_jump= length(x);
for n=1:10
    r=round(50*rand(1))+20;        %% between 20-70
    k=mod(n,2);
    if k==0                        %% k even
        x=[x,ones(1,r)];
    else                           %% k odd
        x=[x,-ones(1,r)];
    end
    index_jump=[index_jump length(x)];
end

L=length(x);                       %% x is a sequence from 0 and 1
t=0:L-1;
y=sin(2*pi/100*t);                 %% create a sinuosly wave with length of L
z=0.15*x + y;                      %% sum .15 or -.15 with sinosal wave (number of discontuonous's point is 10)

Nfft=1024;
step=2*pi/Nfft;
t=-pi:step:pi-step;

%%%% fft of signals 
figure;
subplot(311)
Ffft_rect=abs(fftshift(fft(.15*x,Nfft))); 
plot(t,Ffft_rect); title('fft of .15*x');
subplot(312)
Ffft_sin=abs(fftshift(fft(y,Nfft)));
plot(t,Ffft_sin); title('fft of sin');
subplot(313)
Ffft_sin=abs(fftshift(fft(z,Nfft)));
plot(t,Ffft_sin); title('fft of  z=sin + 0.15x');


%%%%%%%%%%%%%%%
figure;
plot(.15*x);
figure;
plot(y);
figure;
plot(z,'LineWidth',2);
axis([0 length(z) -2 2]);

%%%%%%%%%%%%%  B  %%%%%%%%%%%%%%%
%%%%%%%% a)

%%%%%%%%%%%%%%     Decomposition    %%%%%%%%%%%%%%
for i=0:3
    
    str = sprintf('db%d',2^i);
    
    [c,l]=wavedec(z,4, str );
    
    d{4*i+1}=c(l(1)+l(2)+l(3)+l(4)+1:l(1)+l(2)+l(3)+l(4)+l(5));
    d{4*i+2}=c(l(1)+l(2)+l(3)+1:l(1)+l(2)+l(3)+l(4));
    d{4*i+3}=c(l(1)+l(2)+1:l(1)+l(2)+l(3));
    d{4*i+4}=c(l(1)+1:l(1)+l(2));
    a_4{i+1}=c(1:l(1));
    
end
%%%%%%%%%%%%

%%% B b)
%%%%% plot d1 for 4 Decomposition
q=1;
figure;
for i=0:3
    temp=d{4*i+q};
    subplot(4,1,i+1);
    stem(temp);
    str = sprintf('d-%d - db%d',q,2^i);
    title(str)
end

%%% calculate jump index
temp=d{13};  %%% d1 db8
r= find(temp > .01);
r=r(r>5);
r=[r(1) r];
k=1;
No=0;
temp2=[0];
for i=1:length(r)-1
    
    temp1=r(i+1)-r(i);
    if temp1 > 5
        temp2(k)= r(i);
        k=k+1;
    end
end
index=[temp2 r(end)]
index=[2*temp2 2*r(end)]
index_enhancement=index-6
index_jump
%%%%%%%%%%%%%

%%%% 2 b)

%%%% plot a coefficents for 4 Decomposition
figure;
for i=1:4
    temp=a_4{i};
    subplot(4,1,i);
    plot(temp);
    str = sprintf('a-4 - db%d',2^(i-1));
    title(str)
end



%%%%% C
%%%% plot fft


%%% u=0 => db1 ** u=1 => db2 ** u=2 => db4 ** u=3 => db8
for u=0:3
    q=4*u;
    
    %%% a=1 => d_1 **  a=2 => d_2 **  a=3 => d_3 **  a=4 => d_4 **  a=0 => a_4
    for a=0:4
        
        
        switch a
            case 1 %% d1
                d_1_New=zeros(1,length(d{q+1}));
                
                C=[a_4{u+1}  d{q+4}  d{q+3}  d{q+2} d_1_New];
                L=[length(a_4{u+1})  length(d{q+4})  length(d{q+3})  length(d{q+2}) length(d_1_New) length(z)];
                db_type = sprintf('db%d',2^u);
                z_New=waverec(C,L, db_type);
                
                title_data=sprintf('%s & removed => d-%d',db_type,a);
            case 2 %% d2
                d_2_New=zeros(1,length(d{q+2}));
                
                C=[a_4{u+1}  d{q+4}  d{q+3}  d_2_New d{q+1}] ;
                L=[length(a_4{u+1})  length(d{q+4})  length(d{q+3}) length(d_2_New) length(d{q+1})  length(z)];
                db_type = sprintf('db%d',2^u);
                z_New=waverec(C,L, db_type);
                
                title_data=sprintf('%s & removed => d-%d',db_type,a);
            case 3 %%d3
                d_3_New=zeros(1,length(d{q+3}));
                
                C=[a_4{u+1}  d{q+4}  d_3_New d{q+2}  d{q+1} ];
                L=[length(a_4{u+1})  length(d{q+4}) length(d_3_New) length(d{q+2})  length(d{q+1})  length(z)];
                db_type = sprintf('db%d',2^u);
                z_New=waverec(C,L, db_type);
                
                title_data=sprintf('%s & removed => d-%d',db_type,a);
            case 4 %% d4
                d_4_New=zeros(1,length(d{q+4}));
                
                C=[a_4{u+1}  d_4_New d{q+3}  d{q+2}  d{q+1}];
                L=[length(a_4{u+1})  length(d_4_New) length(d{q+3})  length(d{q+2})  length(d{q+1})  length(z)];
                db_type = sprintf('db%d',2^u);
                z_New=waverec(C,L, db_type);
                
                title_data=sprintf('%s & removed => d-%d',db_type,a);
            case 0 %% c4
                a_4_New=zeros(1,length(a_4{u+1}));
                
                C=[a_4_New d{q+4} d{q+3}  d{q+2}  d{q+1}];
                L=[length(a_4_New)  length(d{q+4}) length(d{q+3})  length(d{q+2})  length(d{q+1})  length(z)];
                db_type = sprintf('db%d',2^u);
                z_New=waverec(C,L, db_type);
                
                title_data=sprintf('%s & removed => a_4',db_type);
        end
        
        figure('units','normalized','outerposition',[0 0 1 1])
        subplot(411)
        plot(z);  title('z (sinuously + rectangular pulse)')
        subplot(412)
        plot(z_New);  str=sprintf('reconstructed signal from z (%s)',title_data);  title(str);
        subplot(413);
        plot(y);      title('sinuously')
        subplot(414);
        plot(.15*x);      title('random rectangular pulse')
        
    end
end




%%%%%% extra code
p=2;
while p==1
    %%% fft plot
    % Nfft=1024;
    % step=2*pi/Nfft;
    % t=-pi:step:pi-step;
    %
    % figure;
    % subplot(211);
    % plot(t, abs(fftshift(fft( y ,Nfft))) );
    % title('sinuosly Signal')
    %
    %
    %
    % subplot(212);
    % plot(t, abs(fftshift(fft( z ,Nfft))) );
    % title('Sinuosly + Rectangular Pulse ')
    
    
    % for i=1:20
    %     temp=[];
    %     r = mod(i,5);
    %     s = floor(i/5);
    %     f=i-s;
    %
    %
    %     if r == 1
    %         figure;
    %     end
    %
    %     if r == 0
    %         subplot(5,1,5);
    %         Ffft_a4{s}=abs(fftshift(fft(a_4{s},Nfft)));
    %         temp=Ffft_a4{s};
    %         str = sprintf('c-4 - db%d',2^(s-1));
    %     else
    %         subplot(5,1,r);
    %         Ffft_d{f}=abs(fftshift(fft(d{f},Nfft)));
    %         temp=Ffft_d{f};
    %         str = sprintf('d-%d - db%d',r,2^s);
    %     end
    %
    %     plot(t,temp);
    %     title(str)
    %     p=2;
    % end
    
    %%%% recounstration and eliminate rectunguale pulse
    


%%%%% decompose with other wavelet

% db_type='db16';
% [c,l]=wavedec(z,4, db_type );
% 
% d_other{1}=c(l(1)+l(2)+l(3)+l(4)+1:l(1)+l(2)+l(3)+l(4)+l(5));
% d_other{2}=c(l(1)+l(2)+l(3)+1:l(1)+l(2)+l(3)+l(4));
% d_other{3}=c(l(1)+l(2)+1:l(1)+l(2)+l(3));
% d_other{4}=c(l(1)+1:l(1)+l(2));
% a_4_other{1}=c(1:l(1));
% 
% for a=0:4
%     
%         switch a
%             case 1 %% d1
%                 d_1_New=zeros(1,length(d_other{1}));
%                 
%                 C=[a_4_other{1}  d_other{4}  d_other{3}  d_other{2} d_1_New];
%                 L=[length(a_4_other{1})  length(d_other{4})  length(d_other{3})  length(d_other{2}) length(d_1_New) length(z)];
%                 
%                 z_New=waverec(C,L, db_type);
%                 
%                 title_data=sprintf('%s & removed => d-%d',db_type,a);
%             case 2 %% d2
%                 d_2_New=zeros(1,length(d{2}));
%                 
%                 C=[a_4_other{1}  d_other{4}  d_other{3}  d_2_New d_other{1}] ;
%                 L=[length(a_4_other{1})  length(d_other{4})  length(d_other{3}) length(d_2_New) length(d_other{1})  length(z)];
%                 
%                 z_New=waverec(C,L, db_type);
%                 
%                 title_data=sprintf('%s & removed => d-%d',db_type,a);
%             case 3 %%d3
%                 d_3_New=zeros(1,length(d{3}));
%                 
%                 C=[a_4_other{1}  d_other{4}  d_3_New d_other{2}  d_other{1} ];
%                 L=[length(a_4_other{1})  length(d_other{4}) length(d_3_New) length(d_other{2})  length(d_other{1})  length(z)];
%              
%                 z_New=waverec(C,L, db_type);
%                 
%                 title_data=sprintf('%s & removed => d-%d',db_type,a);
%             case 4 %% d4
%                 d_4_New=zeros(1,length(d{4}));
%                 
%                 C=[a_4_other{1}  d_4_New d_other{3}  d_other{2}  d_other{1}];
%                 L=[length(a_4_other{1})  length(d_4_New) length(d_other{3})  length(d_other{2})  length(d_other{1})  length(z)];
%                 
%                 z_New=waverec(C,L, db_type);
%                 
%                 title_data=sprintf('%s & removed => d-%d',db_type,a);
%             case 0 %% c4
%                 a_4_New=zeros(1,length(a_4_other{1}));
%                 
%                 C=[a_4_New d_other{4} d_other{3}  d_other{2}  d_other{1}];
%                 L=[length(a_4_New)  length(d_other{4}) length(d_other{3})  length(d_other{2})  length(d_other{1})  length(z)];
%                 
%                 z_New=waverec(C,L, db_type);
%                 
%                 title_data=sprintf('%s & removed => a_4',db_type);
%         end
%         
%         figure('units','normalized','outerposition',[0 0 1 1])
%         subplot(411)
%         plot(z);  title('z (sinuously + rectangular pulse)')
%         subplot(412)
%         plot(z_New);  str=sprintf('reconstructed signal from z (%s)',title_data);  title(str);
%         subplot(413);
%         plot(y);      title('sinuously')
%         subplot(414);
%         plot(.15*x);      title('random rectangular pulse')
%         
% end

end




