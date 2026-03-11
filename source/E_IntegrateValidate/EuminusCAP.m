function good_r = EuminusCAP(params,mflds,pulse4D,Eu)
% CAP for the chebyshev coefficients of E^u_- found in chebInt

%We may work with lower order here (but still keep some tail of 0s for proof)
phi_cheb = pulse4D.a1(1:Eu.nonzero)'; 

params.rho_error = get_rho(params,mflds,pulse4D,phi_cheb,Eu.nonzero);

Lbvp = pulse4D.Lbvp;

ord = params.Eu.order;
nz_ord = Eu.nonzero;
phi_chebN = [phi_cheb;zeros(ord - nz_ord,1)];


a = zeros(4*ord,1);
ICvec = zeros(4,1);

for i = 1:4 
    a((i-1)*ord+1:(i-1)*ord + nz_ord) = Eu.U_1_cheb(:,i);%Defining \bar{a}
    ICvec(i) = chebSum(Eu.U_1_cheb(:,i),-1);%Defining the initial condition vector
end

Ad_N = chebDF(phi_chebN,ord,params,Lbvp,nz_ord);
A_N = Ad_N^-1;
a_bar = intval(1)*a;


Y0 = computeY0(A_N,a_bar,phi_chebN,params,ord,ICvec,params.del,nz_ord,Lbvp);
Y0hat = computeY0hat(A_N,params.rho_error,params.Lbvp,params.del,params.nu,a_bar(1:ord),phi_chebN);
fprintf('Y bounds computed, Y0 = %d, Y0hat = %d\n', sup(Y0), sup(Y0hat));

Z0 = computeZ0(A_N,Ad_N,ord,params.del);
Z1 = computeZ1(A_N,ord,phi_chebN,params.del,params);
Z2hat = computeZ2hat(A_N,pulse4D.Lbvp,params.del,params.nu,params.rho_error);
fprintf('Z bounds computed, Z0 = %d, Z1 = %d, Z2hat = %d\n', mid(Z0), mid(Z1), mid(Z2hat));

good_r = sup((Y0 + Y0hat)/(1-Z0-Z1-Z2hat)); %Radii poly is linear

end