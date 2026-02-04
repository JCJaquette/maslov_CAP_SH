function [params, phi, mani_coeffs] = getparamsBefore(n)

if n == 1

    params.mu=0.05;
    params.nu=1.6;
    params.new_L = 61.421248010101590-3.37;
    params.scale = 3e-1;
    load('test/test_ValidatePulses/mflds1.mat')
    params.mfld_error = 3.2e-11;
    load('test/test_IntegrateValidate/validatedpulse1.mat')


elseif n == 2

    params.mu=0.05;
    params.nu=1.6;
    params.new_L = 61.421248010101590-5.29;
    params.scale = 2.5e-1;
    load('test/test_ValidatePulses/mflds2.mat')
    params.mfld_error = 7.7e-13;
    load('test/test_IntegrateValidate/validatedpulse2.mat')

elseif n ==3

    params.mu=0.2;
    params.nu=1.6;
    params.new_L = 26.181640966137273-11.69;
    params.scale = 3e-1;
    load('test/test_ValidatePulses/mflds3.mat')
    params.mfld_error = 1.9e-18;
    load('test/test_IntegrateValidate/validatedpulse3.mat')

end

mani_coeffs = mflds.unstable.coeffs; %left in natural coordinates

phi = [new_y.phi1, new_y.phi2];

[~, values]= getJacEigs(0, params);
params.lambda = values.u;

end

