function out = taylorSum2D(coeffs,x1,x2)
% taylor coefficients(size order x order) and 2D point in, series at that point out
    
    order = length(coeffs(1,:))-1;
    N = [];
    for i = 0:order    
        N = [N;0:order];       
    end
    M = N';

    x1x2Mat = (x1.^M).*(x2.^N);
    series = x1x2Mat.*coeffs;

    out = sum(sum(series));

end

