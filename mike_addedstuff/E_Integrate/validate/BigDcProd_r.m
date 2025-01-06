function [outmat] = BigDcProd_r(a)
% gives bigger jacobian of a*h wrt h

    N = length(a);

    size = 3*N;

    outmat = zeros(size);

    for i = N+1:size
        b = 0*a;
        b(i) = 1;
        for j = 1:N

            outmat(j,i) = chebcProd_nth(a',b',j-1);%there's probably a better way than this

        end
    end

end

