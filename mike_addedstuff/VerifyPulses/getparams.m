function [params,mflds,psoln] = getparams(n)


if n == 1

    params.rho = 1 - .01;
    params.scale = 3e-1;
    params.mu=0.05;
    params.nu=1.6;
    params.lambda = 0;
    params.cheb.order=450;
    params.mfld.order=40;
    params.tol=4e-16;
    params.L = 0;

    params.bd_scale = .05;
    
    load("psoln1.mat");
    psoln = psoln1;

    params.new = 1.05;

    load('mflds1.mat');

    load('mError1');
    mflds.stable.error = mError1;
    mflds.unstable.error = mError1;

    return

elseif n == 2

    params.rho = 1 - .01;
    params.scale = 2.5e-1;
    params.mu=0.05;
    params.nu=1.6;
    params.lambda = 0;
    params.cheb.order=450;
    params.mfld.order=35;
    params.tol=4e-16;
    params.L = 0;

    params.bd_scale = .08;
    
    load("psoln2.mat");
    psoln = psoln2;

    params.new = 1.05;

    load('mflds2.mat');

    load('mError2');
    mflds.stable.error = mError2;
    mflds.unstable.error = mError2;

    return

elseif n == 3

    params.rho = 1 - .01;
    params.scale = 3e-1;
    params.mu=0.2;
    params.nu=1.6;
    params.lambda = 0;
    params.cheb.order=2^10;
    params.mfld.order=25;
    params.tol=4e-15;
    params.L = 0;

    params.bd_scale = .1;
    
    load("psoln3.mat");
    psoln = psoln3;

    params.new = 1.01;

    load('mflds3.mat');

    load('mError3');
    mflds.stable.error = mError3;
    mflds.unstable.error = mError3;

    return

end

    error('Make sure n=1, 2, or 3')

end

