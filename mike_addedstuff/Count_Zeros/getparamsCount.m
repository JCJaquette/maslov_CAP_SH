function [params,h_cheb,phi_cheb,phiPrime_cheb] = getparamsCount(n) %#ok<STOUT>

if n == 1

    params.mu=0.05;
    params.nu=1.6;
    params.L = 3.37;

    load('Count_Zeros/savedThings/varbs1.mat');

elseif n == 2

    params.mu=0.05;
    params.nu=1.6;
    params.L = 5.29;

    load('Count_Zeros/savedThings/varbs2.mat');

elseif n == 3

    params.mu=0.2;
    params.nu=1.6;
    params.L = 11.69;

    load('Count_Zeros/savedThings/varbs3.mat');

else

    error('Make sure n=1, 2, or 3')

end


    h1coeffs = [h_cheb(1,1);2*h_cheb(2:end,1)];
    h2coeffs = [h_cheb(1,2);2*h_cheb(2:end,2)];
    h4coeffs = [h_cheb(1,4);2*h_cheb(2:end,4)];
    phiP1coeffs = [phiPrime_cheb(1,1);2*phiPrime_cheb(2:end,1)];
    phiP2coeffs = [phiPrime_cheb(1,2);2*phiPrime_cheb(2:end,2)];
    phiP3coeffs = [phiPrime_cheb(1,3);2*phiPrime_cheb(2:end,3)];

    params.f_error = E_phi*(norm(h1coeffs - h4coeffs,1)) ...
                    +E_h*(norm(phiP1coeffs - phiP2coeffs,1));

    params.df_error = E_phi*(norm(h1coeffs - h2coeffs,1)) ...
                    +E_h*(norm(phiP1coeffs - phiP3coeffs,1));

end

