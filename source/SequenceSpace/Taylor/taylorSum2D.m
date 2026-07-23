function out = taylorSum2D(coeffs,x1,x2)
% taylor coefficients(size order x order) and 2D point in, series value at that point out
    
    isIntval = isintval(coeffs) || isintval(x1) || isintval(x2);
    order = length(coeffs(1,:))-1;
    N = 0:order;
    M = N';

    x1M = x1.^M;
    if any(isnan(x1M(:)))
        if isIntval
            x1M = intval(zeros(order+1,1));
        else
            x1M = zeros(order+1,1);
        end

        for m = 0:order
            x1M(m+1) = x1^m;
        end
    end
    
    x2N = x2.^N;
    if any(isnan(x2N(:)))
        if isIntval
            x2N = intval(zeros(1,order+1));
        else
            x2N = zeros(1,order+1);
        end

        for n = 0:order
            x2N(n+1) = x2^n;
        end
    end


    x1x2Mat = x1M * x2N;    
    series = x1x2Mat.*coeffs;
    out = sum(sum(series));

end

