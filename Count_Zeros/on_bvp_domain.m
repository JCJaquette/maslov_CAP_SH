% Count zeros of the determinant on [-L_bvp, L_bvp]

close all
clear

global_min = -1;
global_max = 1;
domINT = infsup(global_min,global_max);

[params,h_cheb,phi_cheb,phiPrime_cheb] = getparamsCount(2);

f_error = infsup(-params.f_error,params.f_error);
df_error = infsup(-params.df_error,params.df_error);

A_cheb = chebstar2fft(h_cheb(:,1),phiPrime_cheb(:,4)) ...
    - chebstar2fft(h_cheb(:,4),phiPrime_cheb(:,1));
APrime_cheb = chebstar2fft(h_cheb(:,1),phiPrime_cheb(:,2)) ...
    - chebstar2fft(h_cheb(:,2),phiPrime_cheb(:,1));

f = @(x) chebSum(A_cheb,x);
df = @(x) chebSum(APrime_cheb,x);

tol = 1e-5;

[count0s, flag] = getZeroCount(domINT, f, f_error, df, df_error, tol);

count0s