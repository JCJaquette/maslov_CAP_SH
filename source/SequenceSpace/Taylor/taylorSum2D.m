function out = taylorSum2D(coeffs,x1,x2)
% taylor coefficients(size order x order) and 2D point in, series value at that point out
%
% Evaluates  out = sum_{q,n} coeffs(q+1,n+1) * x1^q * x2^n.
%
% CAREFUL: do not write  x1.^N'  to build the column of powers.  In MATLAB
% .^ and ' share precedence and associate left-to-right, so x1.^N' parses as
% (x1.^N)' -- a *conjugate* transpose.  With x2 = conj(x1) that silently
% replaces x1^q by conj(x1)^q = x2^q, i.e. it evaluates the wrong series
% (and returns a complex value for a real manifold).  Use an explicit column
% of exponents instead.
    zero_type = 0*x1;

    order = length(coeffs(1,:))-1;
    Nrow = 0:order;         % row of exponents,    for x2
    Ncol = (0:order).';     % column of exponents, for x1

    % The scalar loop is only needed when an interval bracket contains 0,
    % where the vectorized 0^0 would produce NaN.
    if isintval(x1) && in(0,x1)
        x1M = zeros(order+1,1)*zero_type; %convert to intval if necessary
        for m = 0:order
            x1M(m+1) = x1^m;
        end
    else
        x1M = x1.^Ncol;
    end

    if isintval(x2) && in(0,x2)
        x2N = zeros(1,order+1)*zero_type; %convert to intval if necessary
        for n = 0:order
            x2N(n+1) = x2^n;
        end
    else
        x2N = x2.^Nrow;
    end

    x1x2Mat = x1M * x2N;
    series = x1x2Mat.*coeffs;
    out = sum(sum(series));

end
