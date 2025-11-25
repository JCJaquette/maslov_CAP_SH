function [vec1,vec2] = getEu_minusL(mfld_coeffs,phi)
%compute the manifold derivatives wrt sigma1 and sigma2

    d1mfld_coeffs = mfld_coeffs;
    d2mfld_coeffs = mfld_coeffs;
    ord = length(mfld_coeffs(:,1,1));

    for i = 1:ord

        d1mfld_coeffs(i,:,:) = (i-1)*d1mfld_coeffs(i,:,:);
        d2mfld_coeffs(:,i,:) = (i-1)*d2mfld_coeffs(:,i,:);

    end

    d1mfld_coeffs = d1mfld_coeffs(2:end,:,:);
    d2mfld_coeffs = d2mfld_coeffs(:,2:end,:);

    vec1 = get_manifold_point(d1mfld_coeffs,phi(1),phi(2),24);
    vec2 = get_manifold_point(d2mfld_coeffs,phi(1),phi(2),24);

end

