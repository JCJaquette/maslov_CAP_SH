function [c] = make_c(a1,a2,a3,a4,b,params,N)
% makes c from thesis p231

    nz = params.nonzero;

    ba1 = chebstar2fft(b,a1);
    ba1 = [ba1(1:2*nz); zeros(N - 2*nz, 1)];
    bba1 = chebstar2fft(b,ba1);
    bba1 = [bba1(1:3*nz); zeros(N - 3*nz, 1)];
    
    if params.isIntval
        c = intval(zeros(4,N));
    else
        c = zeros(4,N);
    end



    c(1,:) = a4;
    c(2,:) = a3 - 2*a4;
    c(3,:) = -(1+params.mu)*a1 + 2*params.nu*ba1(1:N) - 3*bba1(1:N);
    c(4,:) = a2;


end