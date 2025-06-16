function Y0 = computeY0(A,a_bar,b,params,N,v,del)

    finitepart = A*chebF_intval(a_bar,v,b,N,params);

    Y0s = zeros(1,4)*intval(0);

    for i = 1:4
        Y0s(1) = vectorDelta1norm(finitepart(1:N),del);
        Y0s(2) = vectorDelta1norm(finitepart(N+1:2*N),del);
        Y0s(3) = vectorDelta1norm(finitepart(2*N+1:3*N),del);
        Y0s(4) = vectorDelta1norm(finitepart(3*N+1:4*N),del);
    end

    a_bar3 = a_bar(2*N+1:3*N);
    prod_2 = chebstar2fft_intval(b,a_bar3);
    prod_2 = prod_2(1:3*N);
    %prod_2 = chebstar2(b,a_bar3,3*N);
    prod_3new = chebstar2fft_intval(b,prod_2);
    prod_3new = prod_3new(1:3*N);
    %prod_3new = chebstar2(b,prod_2,3*N); uncomment non-fft version for
    %pulse3/better bounds
    y3tail = params.nu * prod_2 - 3*prod_3new;
    
    for i = 1:3*N

        if i < N
            y3tail(i) = 0;
        else
            y3tail(i) = y3tail(i)/(i+1);
        end

    end

    y3tail_sum = vectorDelta1norm(y3tail,del);

    Y0s(3) = Y0s(3) + y3tail_sum;

    Y0 = max(Y0s);

end

