function [Vs_coeff,Vu_coeff] = getLinSolnCoords(params,bndl,Eucheb,pulse4D)

S = [1, 0, 0, 0; 
     0, 0, 1, 0;
     0, 2, 0, 1;
     0, 1, 0, 0];
sig0 = get_sig0(params,pulse4D);
W_sig0 = bndl_one_point(sig0(1),sig0(2),bndl,params);

U_L = zeros(4,1);
for i = 1:4
    U_L(i) = chebSum(Eucheb(:,i),1);
end

x = W_sig0\(S\U_L);

Vs_coeff = x(1:2);
Vu_coeff = x(3:4);

end