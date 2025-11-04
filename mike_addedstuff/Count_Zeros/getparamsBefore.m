function [params, phi, mani_coeffs] = getparamsBefore(n)

if n == 1

    params.mu=0.05;
    params.nu=1.6;
    params.new_L = 61.421248010101590-3.37;%see L_minus in results folder
    params.scale = 3e-1;
    load('ValidatePulses/saved_things/mflds1.mat')
    params.mfld_error = 3.2e-11;
    load('E_Integrate/pulses/verifiedpulse1.mat')


elseif n == 2

    params.mu=0.05;
    params.nu=1.6;
    params.new_L = 61.421248010101590-5.29;
    params.scale = 2.5e-1;
    load('ValidatePulses/saved_things/mflds2.mat')
    params.mfld_error = 7.7e-13;
    load('E_Integrate/pulses/verifiedpulse2.mat')

elseif n ==3

    params.mu=0.2;
    params.nu=1.6;
    params.new_L = 26.181640966137273-11.69;
    params.scale = 3e-1;
    load('ValidatePulses/saved_things/mflds3.mat')
    params.mfld_error = 1.9e-18;
    load('E_Integrate/pulses/verifiedpulse3.mat')

end

mani_coeffs = mflds.unstable.coeffs; %left in natural coordinates

phi = [new_y.phi1, new_y.phi2];

[~, values]= getJacEigs(0, params);
params.lambda = values.u;

end

