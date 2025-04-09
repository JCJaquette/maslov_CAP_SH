function out = computeZ2hat(A,L,del,nu,Rho)

    N = max(size(A))/4;
    Cnorm = L * (1 + del^2) * 8*(nu*Rho + 6*Rho^2);

    for i = 4:-1:1
        norms(i) = matrixDelta1norm_intval(A((i-1)*N+1:N,2*N+1:3*N),del);
    end

    out = Cnorm * max(norms);

end

