function out = computeZ2hat(A,L,del,nu,Rho)
%See lemma 8.9

    N = max(size(A))/4;
    Cnorm = L * (1 + del^2) * 8*(nu*Rho + 6*Rho^2);%Defined in/below eqn 8.12

    for i = 4:-1:1
        norms(i) = matrixDelta1norm(A((i-1)*N+1:N,2*N+1:3*N),del);
    end

    out = Cnorm * max(norms);

end

