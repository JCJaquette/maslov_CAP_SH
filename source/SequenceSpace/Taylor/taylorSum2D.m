function out = taylorSum2D(coeffs,x1,x2)
% taylor coefficients(size order x order) and 2D point in, series value at that point out
    
    order = length(coeffs(1,:))-1;
    N = 0:order;
    M = N';

    try
        x1M = x1.^M;
    catch
        x1M = zeros(order+1,1);
        for m = 0:order
            x1M(m+1) = x1^m;
        end
    end
    
    try
        x2N = x2.^N;
    catch
        x2N = zeros(1,order+1);
        for n = 0:order
            x2N(n+1) = x2^n;
        end
    end

    x1x2Mat = x1M * x2N;
    series = x1x2Mat.*coeffs;
    out = sum(sum(series));

end

