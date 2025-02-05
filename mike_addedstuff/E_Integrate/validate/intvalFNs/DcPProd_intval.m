function [outmat] = DcPProd_intval(a,b)
% gives jacobian of a*b*h wrt h

    N = length(a);

    c = chebstar2(a,b,N);

    outmat = DcProd_intval(c);

end

