function coeffsOut = hPrimeCoeffs(h_cheb, phi_cheb, params)
% h should be Nx4, function puts out h' coeffs from linear ODE

    phih1_cheb = chebstar2fft(h_cheb(:,1), phi_cheb);
    phi2h1_cheb = chebstar2fft(phi_cheb, phih1_cheb);

    coeffsOut(:,1) = h_cheb(:,4);
    coeffsOut(:,2) = h_cheb(:,3) - 2*h_cheb(:,4);
    coeffsOut(:,3) = 2*params.nu*phih1_cheb - 3*phi2h1_cheb - (1+params.mu)*h_cheb(:,1);
    coeffsOut(:,4) = h_cheb(:,2);

end

