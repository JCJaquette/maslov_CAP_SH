function out = taylorSum2D(coeffs,x1,x2)
% taylor coefficients(size order x order) and 2D point in, series value at that point out
    zero_type = 0*x1;

    order = length(coeffs(1,:))-1;
    N = 0:order;

    if isintval(x1) && ~in(0,x1) % not interval, or x1 doesn't contain 0
        x1M = x1.^N';
    else
        x1M = zeros(order+1,1)*zero_type; %convert to intval if necessary
        for m = 0:order
            x1M(m+1) = x1^m;
        end
    end
    
    if isintval(x2) && ~in(0,x2) % not interval, or x2 doesn't contain 0
        x2N = x2.^N;
    else
        x2N = zeros(1,order+1)*zero_type; %convert to intval if necessary
        for n = 0:order
            x2N(n+1) = x2^n;
        end
    end

    x1x2Mat = x1M * x2N;
    if anynan(x1x2Mat)
        disp('danger!!!!!!!!')
    end
    series = x1x2Mat.*coeffs;
    out = sum(sum(series));

end

