function out = computeZ1(A,N,b,del,params)

    %to bound \pi_N A_N \pi_N D\psi \pi_\infty 
    %use this for L: b = A(1:450,1351:1800)*lvec;
    %               vectorDelta1norm(b,1.001) for the matrices, separate
    %               the Ls out
    Anorm = matrix4Ellnorm(A,N,del);
    L_rNorm = params.L/del;
    B_rNorm = 2/del^(N+1);
    Z_rNorm = (del^-1 + del)*(1 + params.mu + 8*params.nu*vectorDelta1norm(b,del) + 48*vectorDelta1norm(b,del)^2);

    % remainderNorms = [B_rNorm, 0, 0, L_rNorm;
    %                   0, B_rNorm, L_rNorm, 2*L_rNorm;
    %                   Z_rNorm, 0, B_rNorm, 0;
    %                   0, L_rNorm, 0, B_rNorm];

    AB_contribution = Anorm*B_rNorm;

    % norm of A*(the L part)

    lvec = zeros(N,1);
    lvec(end) = params.L;
    AL_rs = zeros(4,4,N);

    for i = 1:4
        for j = 1:4

            AL_rs(i,j,:) = A((i-1)*N + 1: i*N, (j-1)*N + 1: j*N)*lvec; 

        end
    end

    Lnorms = zeros(4,4); %going through the norms of A*Lr matrix
    % the first column is all zero => these norms stay 0

    % norm function doesn't work with thing(i,j,:), so i do this
    vec = zeros(1,N);

    for i = 1:4 
        vec(:) = AL_rs(i,4,:); 
        Lnorms(i,2) = vectorDelta1norm(vec,del);
    end

    for i = 1:4
        vec(:) = AL_rs(i,2,:); 
        Lnorms(i,3) = vectorDelta1norm(vec,del);
    end

    for i = 1:4
        vec(:) = AL_rs(i,1,:) - 2*AL_rs(i,2,:);
        Lnorms(i,4) = vectorDelta1norm(vec,del);
    end

    AL_contribution = norm(Lnorms,inf);

    % norm of A*(the Z part)
    % Since the Z part is only the (1,3) element in the matrix of operators 
    % we only need the 3rd column of A for this part

    bigA13 = zeros(3*N);
    bigA23 = bigA13;
    bigA33 = bigA13;
    bigA43 = bigA13;

    bigA13(1:N,1:N) = A(1:N,2*N + 1:3*N);
    bigA23(1:N,1:N) = A(N + 1:2*N,2*N + 1:3*N);
    bigA33(1:N,1:N) = A(2*N + 1:3*N,2*N + 1:3*N);

    for i = N+1:3*N
        bigA33(i,i) = 1/(2*i);
    end

    bigA43(1:N,1:N) = A(3*N + 1:4*N,2*N + 1:3*N);

    shftfwd = diag(ones(1,3*N-1),1);
    shftbkwd = diag(ones(1,3*N-1),-1);
    bigDcPP = BigDcPProd_r(b);
    bigDcP = BigDcProd_r(b);

    bigDc = 2*params.nu*bigDcP - 3*bigDcPP;

    Dcmns = shftbkwd*bigDc;
    Dcpls = shftfwd*bigDc;
    bigLOL = params.L*(shftbkwd + shftfwd)*eye(3*N);
    bigLOL(1,:) = zeros(1,3*N);

    bigZ = -params.L*(Dcmns - Dcpls) - (1+params.mu)*bigLOL;

    AZnorms = zeros(1,4);

    AZnorms(1) = matrixDelta1norm(bigA13*bigZ,del);
    AZnorms(2) = matrixDelta1norm(bigA23*bigZ,del);
    AZnorms(3) = matrixDelta1norm(bigA33*bigZ,del);
    AZnorms(4) = matrixDelta1norm(bigA43*bigZ,del);

    AZ_contribution = max(AZnorms);

    z1a = AB_contribution + AL_contribution + AZ_contribution;

    %to bound \pi_\infty [.5K^{-1} D\psi - I]

    longb = [b, zeros(1,2*N)];
    bigDF = Dphi_forZ1(longb,3*N,params);

    z1b = matrix4Ellnorm(bigDF,3*N,del)/(2*N);

    out = z1a + z1b;

end

