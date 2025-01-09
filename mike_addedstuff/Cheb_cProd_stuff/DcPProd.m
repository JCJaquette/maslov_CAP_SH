function [outmat] = DcPProd(a,b)
% gives jacobian of a*b*h wrt h

    mat1 = DcProd(a);
    mat2 = DcProd(b);

    outmat = mat1*mat2;

end

