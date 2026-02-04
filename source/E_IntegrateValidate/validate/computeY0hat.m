function out = computeY0hat(A,Rho,L,del,nu,a1,b)
%See lemma 8.2
    N = max(size(A))/4;

    normDiff = L * (1 + del^2) * vectorDelta1norm(a1,del) * ...
        8*Rho*(nu + 3*vectorDelta1norm(b,del) + 6*Rho) /del;%calculation right above lemma 8.2

    for i = 4:-1:1 %Do the norm in the second part of Yhat bound
        norms2(i) = matrixDelta1norm(A((i-1)*N+1:i*N,2*N+1:3*N),del) * normDiff;%why is this bigger?
    end                                                                     

    norms1 = 0*norms2;
    for j = 0:3 %Do the finite vector norms in the first part of the Yhat bound
        for i = 0:3 
            norms1(j+1) = norms1(j+1) + Rho*vectorDelta1norm(A(j*N+1:(j+1)*N,i*N+1),del);
        end
    end

    norms = norms1 + norms2;
    out = max(norms);

end

