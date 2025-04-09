function out = computeY0hat(A,Rho,L,del,nu,a1,b)

    N = max(size(A))/4;

    normDiff = L * (1 + del^2) * vectorDelta1norm(a1,del) * ...
        8*Rho*(nu + 3*vectorDelta1norm(b,del) + 6*Rho) /del;

    for i = 4:-1:1
        norms(i) = matrixDelta1norm_intval(A((i-1)*N+1:i*N,2*N+1:3*N),del) * normDiff;
    end

    out = normDiff * max(norms);

end

