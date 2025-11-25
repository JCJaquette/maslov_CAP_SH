function [outmat] = DcPProd_intval(a,b)
% gives jacobian of a*b*h wrt h

    N = length(a);

    c = chebstar2fft_intval(a,b);
    c = c(1:N);

    outmat = DcProd_intval(c);

end

