% params must contain mu, nu, scale, mflds.order
function mflds = get_mflds(params)
    % get scaled eigenvectors 

    [vectors, values]= getJacEigs_explicit(params);
    
    mflds.vectors = vectors;
    mflds.values  = values;

    mflds.stable.coeffs=calc_proj_coeff(values.s, vectors.s, params);
    mflds.unstable.coeffs=calc_proj_coeff(values.u, vectors.u, params);
    
    
end
