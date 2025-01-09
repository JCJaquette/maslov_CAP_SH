function [outmat] = BigDcProd_r(a)
% gives bigger jacobian of a*h wrt h

    N = length(a);
    a = [a,zeros(1,2*N)];

    outmat = DcProd(a);

    outmat(:,1:N) = zeros(3*N,N);
    outmat(N+1:end,N+1:end) = zeros(2*N,2*N);

end

