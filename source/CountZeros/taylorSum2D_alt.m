function out = taylorSum2D_alt(coeffs,x1,x2)
% taylor coefficients(size order x order) and 2D point in, series value at that point out
    
    order = length(coeffs(1,:)) - 1;

    out = 0;

    for m = 0:order
        for n = 0:order
            out = out + coeffs(m+1, n+1) * (x1^m) * (x2^n);
        end
    end

end