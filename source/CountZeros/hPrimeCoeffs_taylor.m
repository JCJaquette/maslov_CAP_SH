function coeffsOut = hPrimeCoeffs_taylor(h_taylor, phi_taylor, params)
% h should be NxNx4, function puts out h' coeffs from linear ODE

    N = length(h_taylor(:,1,1));
    phih1_taylor = chebstar2fft(h_taylor(:,:,1), phi_taylor);
    phi2h1_taylor = chebstar2fft(phi_taylor, phih1_taylor);

    coeffsOut(:,:,1) = h_taylor(:,:,4);
    coeffsOut(:,:,2) = h_taylor(:,:,3) - 2*h_taylor(:,:,4);
    coeffsOut(:,:,3) = 2*params.nu*phih1_taylor(1:N,1:N) - 3*phi2h1_taylor(1:N,1:N) ...
                     - (1+params.mu)*h_taylor(:,:,1);
    coeffsOut(:,:,4) = h_taylor(:,:,2);

end

