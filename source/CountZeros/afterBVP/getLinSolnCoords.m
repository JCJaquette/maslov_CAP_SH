function [Vs_coeff,Vu_coeff] = getLinSolnCoords(params,bndl,Eucheb,pulse4D)

S = [1, 0, 0, 0; 
     0, 0, 1, 0;
     0, 2, 0, 1;
     0, 1, 0, 0];
sig0 = get_sig0(params,pulse4D);
sig0_r = real(sig0(1));
sig0_i = imag(sig0(1));
W_sig0 = bndl_one_point(sig0_r,sig0_i,bndl,params);
W_sig0 = W_sig0 + infsup(-bndl.r_min.sup,bndl.r_min.sup);

U_L = zeros(4,1);
for i = 1:4
    U_L(i) = chebSum(Eucheb(:,i),1);
end

x = W_sig0\(S\U_L);

Vs_coeff = x(1:2);
Vu_coeff = x(3:4);

end