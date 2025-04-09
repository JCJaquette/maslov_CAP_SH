function [params,manifold_u] = getparamsInt(n)
%params for E_Integrate

if n == 1

    params.mu=0.05;
    params.nu=1.6;
    params.L = 3.37;
    params.cheb.order = 800;
    params.del = 1.01;
    % cheb series error is 2.4e-8, mani is 3.2e-11
    params.rho = 2.4e-8 + 3.2e-11;
    
    load('VerifyPulses/saved_things/mflds1.mat')
    manifold_u.coeffs = mflds.unstable.coeffs;
    load('phis1.mat')
    manifold_u.pulseIC_phi = phis;

    return

elseif n == 2

    params.mu=0.05;
    params.nu=1.6;
    params.L = 5.29;
    params.cheb.order = 400;
    params.del = 1.01;
    % cheb series error is 3.5e-9, mani is 7.7e-13
    params.rho = 3.7e-9 + 7.7e-13;

    load('VerifyPulses/saved_things/mflds2.mat')
    manifold_u.coeffs = mflds.unstable.coeffs;
    load('phis2.mat')
    manifold_u.pulseIC_phi = phis;

    return

elseif n == 3

    params.mu=0.2;
    params.nu=1.6;
    params.L = 11.69;
    params.cheb.order = 600;
    params.del = 1.02;
    % cheb series error is 3.5e-13, mani is 1.9e-18
    params.rho = 3.5e-13 + 1.9e-18;

    load('VerifyPulses/saved_things/mflds3.mat')
    manifold_u.coeffs = mflds.unstable.coeffs;
    load('phis3.mat')
    manifold_u.pulseIC_phi = phis;

    return

end

    error('Make sure n=1, 2, or 3')

end



