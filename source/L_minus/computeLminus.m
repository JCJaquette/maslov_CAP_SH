function params = computeLminus(params, Lminus_cands) 
% Todo: Modify the code so that it doesn't need to take Lminus_cands.

% Eq (3.7)
    % r = (1+params.mu^2)^(1/2);    
    r = (1+abs(params.mu))^(1/2);    

    theta = atan2(sqrt(params.mu),-1);
    if theta <pi/2 || theta >pi
        disp('ERROR ERROR: Wrong branch of arctan')
    end

    r_tilde_abs = abs(cos(theta/2));

    % C.f. Prop 3.4
    % upperbound_old = (7/(16*sqrt(r)))^(1/2)*(2*sqrt(2)*sqrt((1+2*sqrt(r))*(1/r^2 + 1)))^(-1);
    upperbound = (1+r)*(1-r_tilde_abs)/( 40*r^3);

    % We need tau/(1-tau) < upperbound

    [~, values]= getBinfEigs(params);
    lam1 = values.u(1);
    lam2 = values.u(2); 
    
    shiftB = B_infinity(params) - real(lam1).*diag(ones(4,1));
    [V,~] = eigs(shiftB); 
    
    K = max(abs(V),[],'all')*max(abs(V^(-1)),[],'all');
    
    theta = (lam1 + lam2);

    Abs_Re_mu = abs(real(lam1 ));
     
    % need to address the tail terms 
    norms = zeros(params.mfld.order, params.mfld.order);
    for i = 1:params.mfld.order 
        for j = 1:params.mfld.order
            norms(i,j) = vecnorm(params.unstable.coeffs(i,j,:));
        end
    end
    supQ = vecnorm(vecnorm(norms));
    
    C = supQ*(2*params.nu + 6*supQ);

    tau_const = K*C/Abs_Re_mu ;

    if tau_const < 1 
        if tau_const/(1-tau_const) < upperbound 
            Lminus_local=0;
            disp('Found L_minus.')
            disp(Lminus_local)
            params.Lminus = Lminus_local; 
            return
        end
    end



    % tau < tau_const * e^{ - |Re \mu | L_- }

    % We need tau / (1-tau) < upperbound
    if upperbound< 1
        % we want tau < upperbound/(1+upperbound)

        L_approx = log(  upperbound/(1+upperbound) / tau_const)/(-Abs_Re_mu);
        L_approx = L_approx*1.0001;
    end

    tau = exp(- Abs_Re_mu* L_approx) *tau_const;

    if tau/(1-tau) < upperbound 
        disp('Found L_minus.')
        disp(L_approx)
        params.Lminus = L_approx; 
        return
    end

    
    smallEnough = 0; 


    
    i = 1;
    while smallEnough == 0
        cand = Lminus_cands(i);
        tau = K*C/Abs_Re_mu*exp(-Abs_Re_mu*cand);
        
        if tau/(1-tau) < upperbound && tau < 1
            smallEnough = 1;
            disp('Found L_minus.')
            disp(Lminus_cands(i))
        else
            i = i + 1;
        end
    end
    params.Lminus = Lminus_cands(i); 
end
