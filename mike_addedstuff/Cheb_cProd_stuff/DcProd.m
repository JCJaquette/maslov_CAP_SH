function [outmat] = DcProd(a)
% gives jacobian of a*h wrt h

    N = length(a);

    outmat = zeros(N);

    for i = 1:N
        b = 0*a;
        b(i) = 1;
        for j = 1:N

            outmat(j,i) = chebcProd_nth(a',b',j-1);%there's probably a better way than this

        end
    end

end

