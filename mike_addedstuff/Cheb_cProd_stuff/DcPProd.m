function [outmat] = DcPProd(a)
% gives jacobian of a*h wrt h

    N = length(a);

    outmat = zeros(N);

    for i = 1:N
        for j = 1:N

            outmat(i,j) = dcProdidbj(a,i-1,j);

        end
    end


end

