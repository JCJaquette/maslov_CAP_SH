function [outmat] = BigDcPProd_r(a)
% gives bigger jacobian of a*h wrt h

    N = length(a);

    size = 3*N;

    outmat = zeros(size);

    for i = 1:N
        for j = N+1:size

            outmat(i,j) = dcProdidbj(a,i-1,j);

        end
    end


end

