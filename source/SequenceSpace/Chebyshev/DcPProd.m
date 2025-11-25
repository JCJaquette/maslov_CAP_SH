function [outmat] = DcPProd(a,b)
% gives jacobian of a*b*h wrt h

    N = length(a);

    c = chebstar2(a,b,N);

    outmat = DcProd(c);

end

