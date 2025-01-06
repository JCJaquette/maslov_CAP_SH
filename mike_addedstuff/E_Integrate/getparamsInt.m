function [params,IC] = getparamsInt(n)

if n == 1

    params.mu=0.05;
    params.nu=1.6;
    params.L = 3.37;
    
    load('mflds1.mat')
    IC = mflds.unstable.coeffs(1,2,:);
    return

elseif n == 2

    params.mu=0.05;
    params.nu=1.6;
    params.L = 5.29;
        
    load('/Users/mike/Documents/GitHub/maslov_CAP_SH/mike_addedstuff/VerifyPulses/saved_things/mflds2.mat')
    IC = mflds.unstable.coeffs(1,2,:);

    return

elseif n == 3

    params.mu=0.2;
    params.nu=1.6;
    params.L = 11.69;

    load('/Users/mike/Documents/GitHub/maslov_CAP_SH/mike_addedstuff/VerifyPulses/saved_things/mflds3.mat')
    IC = mflds.unstable.coeffs(1,2,:);

    return

end

    error('Make sure n=1, 2, or 3')

end



