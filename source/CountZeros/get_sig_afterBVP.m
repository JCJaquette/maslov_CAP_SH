function sigma = get_sig_afterBVP(params,pulse4D,x)

Lambda = get_Lambda(params,'s');

theta = [params.rho * exp(pulse4D.psi*1i);
         params.rho * exp(-pulse4D.psi*1i)];

sigma = [exp(x*Lambda(1)) * theta(1);
       exp(x*Lambda(2)) * theta(2)];

end