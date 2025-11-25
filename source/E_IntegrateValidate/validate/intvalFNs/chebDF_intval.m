function [out] = chebDF_intval(b,N,params)

    alt = intval(1)*ones(1,N-1);
    for k = 1:N-1

        if mod(k,2) == 1
            alt(k) = -alt(k);
        end

    end    

    twos = 2*alt;

    lowright = intval(1)*zeros(N-1,N-1);

    for i = 1:N-1
        lowright(i,i) = intval(2*i);
    end

    diagM = [intval(1),twos;
           zeros(N-1,1),lowright];
    % diagM is derivative of psi_i wrt (a_i), appears as
    % the diagonal blocks
 
    subLOL = intval(1)*zeros(N-1,N);
    
    for i = 1:N-1

        subLOL(i,i) = -params.L;
        if i ~= N-1
            subLOL(i,i+2) = params.L;
        end

    end

    LOL = [intval(1)*zeros(1,N);intval(1)*subLOL];
    % LOL matrix comes up in derivative of psi_i wrt (a_j) with i=/=j  
    % if a_j appears in c_i

    shftfwd = intval(1)*diag(ones(1,N-1),1);
    shftbkwd = intval(1)*diag(ones(1,N-1),-1);    
    DcPP = DcPProd_intval(b,b);
    DcP = DcProd_intval(b);

    Dc = 2*params.nu*DcP - 3*DcPP;
    % derivative of cProds in c_3 wrt a_1

    Dcmns = shftbkwd*Dc;
    Dcpls = shftfwd*Dc;

    C = -params.L*(Dcmns - Dcpls) - (1+params.mu)*LOL;
    C(1,1:N) = zeros(1,N);
    % C is derivative of psi_3 wrt a_1 which is a bit more complicated due
    % to cProds

    zps = zeros(N,N);

    out = [diagM, zps, zps, LOL;
        zps, diagM, LOL, -2*LOL;
        C, zps, diagM, zps;
        zps, LOL, zps, diagM];
    % full DF matrix

end