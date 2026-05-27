function [count0s,flag] = countAfterBVP(params,bndl,mflds,pulse4D,Eu,tol,BOOL_plot)
% Count zeros of the determinant on [-L_conj, -L_bvp]

global_min = 0;
global_max = sup(mflds.Lplus - pulse4D.Lbvp);
domINT = infsup(global_min,global_max);

%get eta,beta,gamma

[tbeta,tgamma] = getLinSolnCoords(params,bndl,Eu.U_1_cheb,pulse4D);

teta = getLinSolnCoords(params,bndl,Eu.U_vpp_cheb,pulse4D);

tildes = [tbeta,tgamma,teta];

f = @(x) get_function(params,bndl,mflds,pulse4D,tildes,x);
df = @(x) 1;
f_error = 0;
df_error = 0;

[count0s, flag] = getZeroCount(domINT, f, f_error, df, df_error, tol, BOOL_plot);

end