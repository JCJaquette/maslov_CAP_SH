function sigma = get_sig_afterBVP(params,pulse4D,x)

Lambda = get_Lambda(params,'s');

psi_int = infsup(pulse4D.psi - pulse4D.r,pulse4D.psi + pulse4D.r);

theta = [params.rho * exp(psi_int*1i);
         params.rho * exp(-psi_int*1i)];

sigma = [exp(x*Lambda(1)) * theta(1);
       exp(x*Lambda(2)) * theta(2)];

end