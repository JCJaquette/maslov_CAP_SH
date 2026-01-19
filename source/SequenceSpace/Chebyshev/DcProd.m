function [outmat] = DcProd(a)
% gives jacobian of a*h wrt h

if isintval(a)
    zero = intval(0);
else
    zero = 0;
end

    N = length(a);

    mat1 = toeplitz(a);

    mat2 = zero*zeros(N);

    for i = 2:N

        mat2(1:N-i+1,i) = a(i:N);

    end

    outmat = mat1 + mat2;

    
end

