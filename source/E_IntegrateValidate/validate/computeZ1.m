function out = computeZ1(A,N,b,del,Lbvp,params)
%See lemmas 8.4-8.8

    %to bound A_N \pi_N D\psi \pi_\infty 

    Anorm = matrix4Ellnorm(A,N,del);

    % norm of A*(the B part)

    B_rNorm = 2/del^(N+1);
    ZB_contribution = Anorm*B_rNorm;

    % norm of A*(the L part)
    
    AL_rs = intval(0)*zeros(4,4,N);
    lvec = zeros(N,1);
    lvec(end) = Lbvp;

    for i = 1:4
        for j = 1:4

            AL_rs(i,j,:) = A((i-1)*N + 1: i*N, (j-1)*N + 1: j*N)*lvec;

        end
    end

    Lnorms = intval(0)*zeros(4,4); %going through the norms of A*Lr matrix
    % the first column is all zero => these norms stay 0

    % norm function doesn't work with thing(i,j,:), so i do this
    vec = intval(0)*zeros(1,N);

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

    ZL_contribution = norm(Lnorms,inf)/del^(N+1);

    % norm of A*(the C part)
    % Since the C part is only the (1,3) element in the matrix of operators 
    % we only need the 3rd column of A for this part

    bigA13 = intval(0)*zeros(3*N);
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

    longb = [b; zeros(2*N,1)];

    bigDcPP = DcPProd(longb,longb);
    bigDcP = DcProd(longb);        %Get the full C operator

    bigDc = 2*params.nu*bigDcP - 3*bigDcPP;
    bigDc(:,1:N) = zeros(3*N,N);             % Since we do A*pi_N C pi_\infty 
    bigDc(N+1:end,N+1:end) = zeros(2*N,2*N); % See lemma 8.6

    Dcmns = shftbkwd*bigDc;
    Dcpls = shftfwd*bigDc;

    bigC = -Lbvp*(Dcmns - Dcpls);
    bigC(N,N+1) = bigC(N,N+1) - Lbvp*(1+params.mu);
    bigC(N+1,N+2) = bigC(N+1,N+2) - Lbvp*(1+params.mu);

    ACnorms = intval(0)*zeros(1,4);

    ACnorms(1) = matrixDelta1norm(bigA13*bigC,del);
    ACnorms(2) = matrixDelta1norm(bigA23*bigC,del);
    ACnorms(3) = matrixDelta1norm(bigA33*bigC,del);
    ACnorms(4) = matrixDelta1norm(bigA43*bigC,del);

    ZC_contribution = max(ACnorms);

    z1a = ZB_contribution + ZL_contribution + ZC_contribution;

    %to bound \pi_\infty [.5K^{-1} D\psi - I]

    bigDF = Dphi_forZ1(longb,3*N,Lbvp,params);

    ZD = matrix4Ellnorm(bigDF,3*N,del)/(2*N);

    out = z1a + ZD;

end

