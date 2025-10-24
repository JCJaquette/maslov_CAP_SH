% params must contain mu, nu, scale, mflds.order
function mflds = get_mflds(params)
    % get scaled eigenvectors 

    [vectors, values]= getJacEigs_explicit(params);
    
    mflds.vectors = vectors;
    mflds.values  = values;

    % Compute the manifolds, with the eigenvectors scaled as desired
    mflds.stable.coeffs=calc_proj_coeff(values.s, vectors.s * params.scale, params);
    mflds.unstable.coeffs=calc_proj_coeff(values.u, vectors.u * params.scale, params);
    
    % Define the validation radius (Default is NaN)
    mflds.stable.r_min = NaN;
    mflds.unstable.r_min = NaN;
    
end
