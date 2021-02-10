function conv_Result = calc_Conv(f,g,t)

l_f=length(f);
% conv_Result=zeros(l_f);
for i=0:16000-1
    temp_f=zeros(1,l_f);
    temp_f(i+1:end)=f(1:end-i);
    
    
    
    temp=temp_f.*g;
    conv_Result(i+1)=trapz(t,temp);
%     a=find(temp);
end

end