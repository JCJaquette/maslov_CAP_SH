function [outmat] = DcProd_intval(a)
% gives jacobian of a*h wrt h

    N = length(a);

    mat1 = toeplitz(a);

    mat2 = intval(1)*zeros(N);

    for i = 2:N

        mat2(1:N-i+1,i) = a(i:N);

    end

    outmat = mat1 + mat2;

    
end
