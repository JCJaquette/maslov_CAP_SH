function [outmat] = BigDcPProd_r_intval(a,b)
% gives bigger jacobian of a*a*h wrt h

    N = length(a);

    a = [a,zeros(1,2*N)];
    b = [b,zeros(1,2*N)];

    outmat = DcPProd_intval(a,b);

    outmat(:,1:N) = zeros(3*N,N);
    outmat(N+1:end,N+1:end) = zeros(2*N,2*N);    

end

