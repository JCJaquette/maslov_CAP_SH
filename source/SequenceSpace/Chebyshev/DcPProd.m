function [outmat] = DcPProd(a,b)
% gives jacobian of a*b*h wrt h

    N = length(a);

    c = chebstar2fft(a,b);
    c = c(1:N);

    outmat = DcProd(c);

end

