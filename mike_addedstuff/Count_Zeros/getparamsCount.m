function [params,h_cheb,phi_cheb,phiPrime_cheb] = getparamsCount(n) %#ok<STOUT>

if n == 1

    params.mu=0.05;
    params.nu=1.6;
    params.L = 3.37;

    load('Count_Zeros/savedThings/varbs1.mat');

    return

elseif n == 2

    params.mu=0.05;
    params.nu=1.6;
    params.L = 5.29;

    load('Count_Zeros/savedThings/varbs2.mat');

    return

elseif n == 3

    params.mu=0.2;
    params.nu=1.6;
    params.L = 11.69;

    load('Count_Zeros/savedThings/varbs3.mat');

    return

end

    error('Make sure n=1, 2, or 3')

end

