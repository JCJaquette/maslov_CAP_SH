function out = matrix4Ellnorm_intval(M,N,del)
% max norm for the (\ell _\delta ^1)^4 space

    norms = intval(0)*zeros(4,4);

    for i = 1:4
        for j = 1:4

            norms(i,j) = matrixDelta1norm_intval(M( (i-1)*N + 1:i*N, (j-1)*N + 1:j*N ),del);

        end
    end

    out = norm(norms,inf);

end