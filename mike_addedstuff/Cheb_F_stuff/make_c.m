function [c] = make_c(a1,a2,a3,a4,b,params,N)
% makes c from thesis p231

    ba1 = 0*a1;
    for i = 1:N
        ba1(i) = chebcProd_nth(a1',b',i-1);
    end

    bba1 = 0*a1;
    for i = 1:N
        bba1(i) = chebcPProd_kth(a1',b',b',i-1);
    end

    c = zeros(4,N);
    c(1,:) = a4;
    c(2,:) = a3 - 2*a4;
    c(3,:) = -(1+params.mu)*a1 + 2*params.nu*ba1 - 3*bba1;
    c(4,:) = a2;


end