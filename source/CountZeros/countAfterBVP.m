function [count0s,flag] = countAfterBVP(params,bndl,mflds,pulse4D,tol,BOOL_plot)
% Count zeros of the determinant on [-L_conj, -L_bvp]

global_min = 0;
global_max = sup(mflds.Lplus - pulse4D.Lbvp);
domINT = infsup(global_min,global_max);


f = @(x) a1(x)*b4(x) - a4(x)*b1(x);

[count0s, flag] = getZeroCount(domINT, f, f_error, df, df_error, tol, BOOL_plot);

end