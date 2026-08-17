function [count_out,flag] = countAfterBVP(params,bndl,mflds,pulse4D,Eu,tol,BOOL_plot)
% Count zeros of the determinant on [-L_conj, -L_bvp]

global_min = 0;
global_max = sup(mflds.Lplus - pulse4D.Lbvp);
domINT = infsup(global_min,global_max);

%get eta,beta,gamma

[tbeta,tgamma] = getLinSolnCoords(params,bndl,Eu.U_1_cheb,pulse4D);

teta = getLinSolnCoords(params,bndl,Eu.U_vpp_cheb,pulse4D);

tildes = [tbeta,tgamma,teta];

sig0 = get_sig0(params,pulse4D);

S = [1, 0, 0, 0; 
     0, 0, 1, 0;
     0, 2, 0, 1;
     0, 1, 0, 0];

f = @(x) get_F_afterBVP(params,bndl,mflds,sig0,pulse4D,tildes,x,S);
df = @(x) get_df_afterBVP(params,bndl,mflds,sig0,pulse4D,tildes,x,S);

[count_out, zero_int, flag] = getZeroCount(domINT, f, 0, df, 0, tol, BOOL_plot);
%note ferror and dferror set to zero since it's all contained in intvals

fprintf('%d zeros found\n', count_out);

if ~isempty(zero_int)

    z_check_l = length(zero_int);

    a1 = @(x) get_a1x(params,bndl,mflds,sig0,pulse4D,tildes,x,S);
    a4 = @(x) get_a4x(params,bndl,mflds,sig0,pulse4D,tildes,x,S);
    b1 = @(x) get_b1x(params,bndl,mflds,sig0,pulse4D,tildes,x,S);
    b4 = @(x) get_b4x(params,bndl,mflds,sig0,pulse4D,tildes,x,S);
    
    for i = 1:z_check_l
    
        fprintf('%s for the %dth zero.\n', ...
    intersectioncheck(zero_int, f, df, a1, a4, b1, b4), i);
    
    end

end

end