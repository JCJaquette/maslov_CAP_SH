function [out] = Dphi_forZ1(b,N,params)
%copy paste of DF but without the diagonal which comes from LHS of the
%functional equation

    alt = ones(1,N-1);
    for k = 1:N-1

        if mod(k,2) == 1
            alt(k) = -1;
        end

    end    

    twos = 2*alt.*ones(1,N-1);

    lowright = zeros(N-1,N-1); %do nothing since we want to take out the K part

    diagM = [1,twos;
           zeros(N-1,1),lowright];
    % diagM is derivative of psi_i wrt (a_i), appears as
    % the diagonal blocks
 
    subLOL = zeros(N-1,N); 
    
    for i = 1:N-1

        subLOL(i,i) = -params.L;
        if i ~= N-1
            subLOL(i,i+2) = params.L;
        end

    end

    LOL = [zeros(1,N);subLOL];
    % LOL matrix comes up in derivative of psi_i wrt (a_j) with i=/=j  
    % if a_j appears in c_i

    shftfwd = diag(ones(1,N-1),1);
    shftbkwd = diag(ones(1,N-1),-1);    
    DcPP = DcPProd(b,b);
    DcP = DcProd(b);

    Dc = 2*params.nu*DcP - 3*DcPP;
    % derivative of cProds in c_3 wrt a_1

    Dcmns = shftbkwd*Dc;
    Dcpls = shftfwd*Dc;

    A = -params.L*(Dcmns - Dcpls) - (1+params.mu)*LOL;
    A(1,1:N) = zeros(1,N);
    % A is derivative of psi_3 wrt a_1 which is a bit more complicated due
    % to cProds

    zps = zeros(N,N);

    out = [diagM, zps, zps, LOL;
        zps, diagM, LOL, -2*LOL;
        A, zps, diagM, zps;
        zps, LOL, zps, diagM];
    % full DF matrix
end

