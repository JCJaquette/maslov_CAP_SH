function APrime_cheb = get_APrime_cheb(a_cheb,b_cheb,phi_cheb,params)
%apply product rule to a1b4 - a4b1

    aprime_cheb = hPrimeCoeffs(a_cheb,phi_cheb,params);
    bprime_cheb = hPrimeCoeffs(b_cheb,phi_cheb,params);

    prodrule1 = chebstar2fft(aprime_cheb(:,1),b_cheb(:,4)) ...
                + chebstar2fft(a_cheb(:,1),bprime_cheb(:,4));
    prodrule2 = chebstar2fft(aprime_cheb(:,4),b_cheb(:,1)) ...
                + chebstar2fft(a_cheb(:,4),bprime_cheb(:,1));

    APrime_cheb = prodrule1 - prodrule2;
    
end

