function [params,manifold_u] = getparamsInt(n)
%params for E_Integrate

if n == 1

    params.mu=0.05;
    params.nu=1.6;
    params.L = 3.37;
    params.cheb.order = 600;
    params.del = 1.01;
    % cheb series error is 2.4e-8, mani is 3.2e-11, max|phi| = 0.996067953847602  
    mani_error = 3.2e-11; ICerror = 2*pi/log(1/0.996067953847601) * mani_error;
    params.rho = max(2.4e-8,ICerror);
    load('ValidatePulses/saved_things/mflds1.mat')
    manifold_u.coeffs = mflds.unstable.coeffs;

    return

elseif n == 2

    params.mu=0.05;
    params.nu=1.6;
    params.L = 5.29;
    params.cheb.order = 600;
    params.del = 1.01;
    % cheb series error is 3.5e-9, mani is 7.7e-13, max|phi| = 0.809718650746043
    mani_error = 7.7e-13; ICerror = 2*pi/log(1/0.809718650746042) * mani_error;  
    params.rho = max(3.5e-9,ICerror);

    load('ValidatePulses/saved_things/mflds2.mat')
    manifold_u.coeffs = mflds.unstable.coeffs;

    return

elseif n == 3

    params.mu=0.2;
    params.nu=1.6;
    params.L = 11.69;
    params.cheb.order = 600;
    params.del = 1.02;
    % cheb series error is 3.5e-13, mani is 1.9e-18, max|phi| = 0.865601439415368
    mani_error = 1.9e-18; ICerror = 2*pi/log(1/0.865601439415367) * mani_error;
    params.rho = max(3.5e-13,ICerror);

    load('ValidatePulses/saved_things/mflds3.mat')
    manifold_u.coeffs = mflds.unstable.coeffs;

    return

end

    error('Make sure n=1, 2, or 3')

end



