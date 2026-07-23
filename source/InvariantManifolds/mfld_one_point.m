function val = mfld_one_point(phi1, phi2, coeff, params)

    x1 = phi1+1i*phi2; x2 = phi1-1i*phi2;
    if isintval(coeff)
        val = intval(1)*zeros(4,1);
    else
        val = zeros(4,1);
    end

    for k = 1:4
        val(k) = taylorSum2D(coeff(:,:,k),x1,x2);
    end

    % if norm(val-real(val)) > 1e-10 
    %      msg = 'Error occurred. The manifold is complex valued.';
    %      error(msg);
    % end %removed because this function is used to compute bundles as well

end 

