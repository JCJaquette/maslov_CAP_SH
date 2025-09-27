function [outputArg1,outputArg2] = bundle_rad_poly(params)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
    mflds.coeff.s = calc_proj_coeff(values.s, vectors.s, params); 
    G_hat = DFQbundle(params, mflds); % I think this is \hat G
end