function L_out = computeLplus(params,bndl,mflds,U_1)

    S = [1, 0, 0, 0; 
         0, 0, 1, 0;
         0, 2, 0, 1;
         0, 1, 0, 0];
    sig0 = mflds.sig0;

    W_sig0 = bndl_one_point(sig0(1),sig0(2),bndl,params);
    %work out error in sig0, bndl
    U1_L = zeros(4,1);
    for i = 1:4
        U1_L(i) = chebSum(U_1(:,i),1);
    end
    
    x = W_sig0\(S\U1_L); %See eqn 6.14 and the one below that

    tbeta = x(1:2); tgamma = x(3:4);

    mu_s = abs(real(mflds.values.u(1)));
    
    V = [mflds.vectors.u, mflds.vectors.s];

    biggestVec = intval(zeros(1,4));
    for i = 1:4
        biggestVec(i) = norm(V(:,1));
    end
    biggestVec = max(biggestVec);

    Vcheck = S*V;

    Vcheck=Vcheck/Vcheck(2,2);  % Renormalize so that 3.6 holds exactly 
    Vnorm = norm(Vcheck);
    Vnorm1 = norm(Vcheck^-1);
    
    for i = 1:4
        component_norms(i)= sum(abs(mflds.unstable.coeffs(:,:,i)),'all');
    end

    p = 1;
    manifold_norm = norm(component_norms,p) + mflds.unstable.r_min;
    
    C = manifold_norm *(2*params.nu + 6*manifold_norm );

    tau = @(L) Vnorm1*Vnorm * C * exp(-mu_s*L)/mu_s; %This bound comes from similar reasoning to eqn 5.1
    %Need tau<1 for a good L, add as condition in theorem
    eps_0 = @(L) tau(L)/(1 - tau(L))*biggestVec; %See 

    eps_beta = @(L) exp(-mu_s*L) * norm(tbeta);
    eps_gamma = @(L) exp(-mu_s*L) / norm(tgamma);

    Vu14 = Vcheck([1,4],1:2);
    Vu12 = Vcheck([1,2],1:2);
    Vs14 = Vcheck([1,4],3:4);
    Vs12 = Vcheck([1,2],3:4);
    Vs34 = Vcheck([3,4],3:4);

    M1 = Vs34'*Vu12 + Vu12'*Vs34;%This is right, fix mistake in paper3
    M2 = inv(Vs14)*Vu14;

    C_M4 = @(L) norm(M2)*(norm(inv(Vu14)) + norm(inv(Vs14)) * (1 + eps(L)*norm(inv(Vu14))) ...
                                                            / (1 - eps(L)*norm(inv(Vs14))));
    C_M3 = @(L) 2*(norm(Vs34) + norm(Vs12)) + 2*eps_0(L);

    MM = (M2'*M1)'*(M2'*M1);
   [v,d] = eig(mid(MM));
   [m_,ind] = min(diag(d));
   [mu,~] = verifyeig(MM,m_,v(:,ind));
   sigmin = sqrt(abs(mu)); %See eqn 3.34 in paper 2


    to_bound = @(L) eps_0(L)*C_M3(L)*norm(M2) + ...
        (eps_0(L)*C_M3(L) + norm(M1))*(eps_0(L)*C_M4(L) + eps_beta(L)*eps_gamma(L));

    L_out = Lplus_bisection(to_bound, @(L) eps_0(L)*norm(inv(Vs14)), sigmin);

end