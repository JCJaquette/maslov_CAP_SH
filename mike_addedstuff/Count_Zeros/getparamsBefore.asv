function [params, phi, lambdas, mani_coeffs] = getparamsBefore(n)

if n == 1

    params.mu=0.05;
    params.nu=1.6;
    params.new_L = 10;
    params.scale = 3e-1;
    load('ValidatePulses/saved_things/mflds1.mat')
    load('E_Integrate/pulses/verifiedpulse1.mat')

elseif n == 2

    params.mu=0.05;
    params.nu=1.6;
    params.new_L = 5.29;
    params.scale = 2.5e-1;
    load('ValidatePulses/saved_things/mflds2.mat')
    load('E_Integrate/pulses/verifiedpulse2.mat')

elseif n ==3

    params.mu=0.2;
    params.nu=1.6;
    params.new_L = 11.69;
    params.scale = 3e-1;
    load('ValidatePulses/saved_things/mflds3.mat')
    load('E_Integrate/pulses/verifiedpulse3.mat')

end

mani_coeffs = mflds.unstable.coeffs;
phi = [new_y.phi1, new_y.phi2];

[~, values]= getJacEigs(0, params);
lambdas = values.u;

end

