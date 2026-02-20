function sig0 = get_sig0(params,pulse4D)
% Get sig0 with error (Here, sig0 are the R^2 coords which get mapped by 
% the stable manifold to \varphi(L^+_bvp), ie \rho e^{\sigma i} below eqn 2.5)

sigma_int = infsup(pulse4D.psi - pulse4D.r,pulse4D.psi + pulse4D.r);

sig0 = [real(params.rho*exp(sigma_int*1i));
        imag(params.rho*exp(sigma_int*1i))];

end