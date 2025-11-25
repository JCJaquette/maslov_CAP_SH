function phiprime_cheb = RHSofODE_coeffs(pulse_cheb,params)
%gives coeffs of eqn 9.4 in thesis

    phiprime_cheb = 0*pulse_cheb;
    ord = length(pulse_cheb);

    % now RHS of SH ODE
    
    phi2 = chebstar2fft(pulse_cheb(1,:),pulse_cheb(1,:));
    phi3 = chebstar2fft(phi2,pulse_cheb(1,:));

    phiprime_cheb(1,:) = pulse_cheb(2,:);
    phiprime_cheb(2,:) = pulse_cheb(3,:);
    phiprime_cheb(3,:) = pulse_cheb(4,:);
    phiprime_cheb(4,:) = -2*pulse_cheb(3,:) - (params.mu + 1)*pulse_cheb(1,:) ...
                        + params.nu*phi2(1:ord) - phi3(1:ord);

    phiprime_cheb = phiprime_cheb*params.L;

end

