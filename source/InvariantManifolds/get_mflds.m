% params must contain mu, nu, scale, mflds.order
function mflds = get_mflds(params)
    % get scaled eigenvectors 
    [vectors, values]= getJacEigs(0, params);
    
    mflds.stable.coeffs=calc_proj_coeff(values.s, vectors.s, params);
    mflds.unstable.coeffs=calc_proj_coeff(values.u, vectors.u, params);
    
    
end
