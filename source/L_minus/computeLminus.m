function L_approx_out = computeLminus(params,mflds) 
    % Stores 

% Eq (3.7) 
    r = (1+abs(params.mu))^(1/2);    

    % C.f. Prop 3.4
    upperbound = 1/( 8*sqrt(r^3));

    % We need tau/(1-tau) < upperbound

    % See (3.14)
    % TODO : Why are these not stored???
    [~, values]= getBinfEigs(params);
    lam1 = values.u(1);
    lam2 = values.u(2); 
    
    shiftB = B_infinity(params);
    [V,~] = eigs(shiftB); 
    
    % (3.14)
    K = max(abs(V),[],'all')*max(abs(V^(-1)),[],'all');
    


    Abs_Re_mu = abs(real(lam1 ));
     
    % TODO need to address the tail terms 
    norms = zeros(params.mfld.order, params.mfld.order);
    for i = 1:params.mfld.order 
        for j = 1:params.mfld.order
            norms(i,j) = vecnorm(mflds.unstable.coeffs(i,j,:)); % What is this??
        end
    end
    supQ = vecnorm(vecnorm(norms));
    
    % (3.3)
    C = supQ*(2*params.nu + 6*supQ);

    % (5.1)
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
        L_approx_out = L_approx;
        params.Lminus = L_approx; 
        return
    else
        L_approx_out = nan;
        disp('Error: Could not find -L_{conj}!')
    end


end
