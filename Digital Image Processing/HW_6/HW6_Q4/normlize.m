function output = normlize(input,type)

[M,N] = size(input);
min_m = min(input(:));
max_m = max(input(:));
d=max_m-min_m;

switch type
    case 1  %%%% for Econt and Ecurv
        if d==0
            output = zeros(M,N);
        else
            output = (input - min_m)./d;
        end
        
    case 2 %%%% for E grad
        if d < 5
            min_m = max_m - 5;
            d=5;
        end
        output = (input - min_m)./ d;
end

end