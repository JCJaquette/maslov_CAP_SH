function [count_out,flag] = countBVP(pulse4D,Euminus,tol,BOOL_plot)
% Count zeros of the determinant on [-L_bvp, L_bvp]

global_min = -1;
global_max = 1;
domINT = infsup(global_min,global_max);

phiPrime_cheb = intval(1)*Euminus.U_vpp_cheb;
h_cheb = intval(1)*Euminus.U_1_cheb;

[f_errbd,df_errbd] = get_det_error(pulse4D,Euminus);

f_error = infsup(-f_errbd,f_errbd);
df_error = infsup(-df_errbd,df_errbd);

A_cheb = chebstar2fft(h_cheb(:,1),phiPrime_cheb(:,4)) ...
    - chebstar2fft(h_cheb(:,4),phiPrime_cheb(:,1));
APrime_cheb = chebstar2fft(h_cheb(:,1),phiPrime_cheb(:,2)) ...
    - chebstar2fft(h_cheb(:,2),phiPrime_cheb(:,1));

f = @(x) chebSum(A_cheb,x);
df = @(x) chebSum(APrime_cheb,x);

[count_out, zero_int, flag] = getZeroCount(domINT, f, f_error, df, df_error, tol, BOOL_plot);

fprintf('%d zeros found\n', count_out);

if ~isempty(zero_int)

    z_check_l = length(zero_int);
    
    a1 = @(t) chebSum(h_cheb(:,1),t);
    a4 = @(t) chebSum(h_cheb(:,4),t);
    b1 = @(t) chebSum(phiPrime_cheb(:,1),t);
    b4 = @(t) chebSum(phiPrime_cheb(:,4),t);
    
    for i = 1:z_check_l
    
        fprintf('%s for the %dth zero.\n', ...
    intersectioncheck(zero_int, f, df, a1, a4, b1, b4), i);
    
    end

end

end