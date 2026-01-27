function good_r = EuminusCAP(params,phi_cheb,U_1_cheb)
% CAP
ord = params.Eu.order;
a = zeros(4*ord,1);
ICvec = zeros(4,1);

for i = 1:4 
    a((i-1)*ord+1:i*ord) = U_1_cheb(:,i);
    ICvec(i) = chebSum(U_1_cheb(:,i),-1);
end

Ad_N = chebDF(phi_cheb,ord,params);
A_N = Ad_N^-1;
a_bar = intval(1)*a;


Y0 = computeY0(A_N,a_bar,phi_cheb,params,ord,ICvec,params.del);
Y0hat = computeY0hat(A_N,params.rho,params.Lbvp,params.del,params.nu,a_bar(1:ord),phi_cheb);
fprintf('Y bounds computed, Y0 = %d, Y0hat = %d\n', sup(Y0), sup(Y0hat));

Z0 = computeZ0(A_N,Ad_N,ord,params.del);
Z1 = computeZ1(A_N,ord,phi_cheb,params.del,params);
Z2hat = computeZ2hat(A_N,params.Lbvp,params.del,params.nu,params.rho);
fprintf('Z bounds computed, Z0 = %d, Z1 = %d, Z2hat = %d\n', mid(Z0), mid(Z1), mid(Z2hat));

rs = 0:10^-12:10^-6;
radii_poly = Y0 + Y0hat - (1-Z0-Z1-Z2hat)*rs;

good_r = sup((Y0 + Y0hat)/(1-Z0-Z1-Z2hat))

end