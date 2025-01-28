function [params,IC] = getparamsInt(n)
%params for E_Integrate

if n == 1

    params.mu=0.05;
    params.nu=1.6;
    params.L = 3.37;
    params.ord = 800;
    params.del = 1.005;
    
    load('VerifyPulses/saved_things/mflds1.mat')
    IC = mflds.unstable.coeffs(1,2,:);
    return

elseif n == 2

    params.mu=0.05;
    params.nu=1.6;
    params.L = 5.29;
    params.ord = 1300;
    params.del = 1.0033;
        
    load('VerifyPulses/saved_things/mflds2.mat')
    IC = mflds.unstable.coeffs(1,2,:);

    return

elseif n == 3

    params.mu=0.2;
    params.nu=1.6;
    params.L = 11.69;
    params.ord = 1700;
    params.del = 1.004;

    load('VerifyPulses/saved_things/mflds3.mat')
    IC = mflds.unstable.coeffs(1,2,:);

    return

end

    error('Make sure n=1, 2, or 3')

end



