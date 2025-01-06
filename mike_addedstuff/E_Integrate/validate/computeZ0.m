function Z0 = computeZ0(A,Ad,N,del)

    B = eye(4*N) - A*Ad;
    Zijs = zeros(4,4);

    for i = 1:4
        for j = 1:4
            Zijs(i,j) = matrixDelta1norm(B((i-1)*N+1 : i*N,(j-1)*N+1 : j*N),del);
        end
    end

    B_rowsums = zeros(4,1);

    for i = 1:4
        B_rowsums(i) = norm(Zijs(i,:),1);
    end

    Z0 = max(B_rowsums);

end

