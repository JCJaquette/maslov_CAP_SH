function [outputArg1,outputArg2] = EuminusCAP(params,y,U_1_cheb)
% CAP
ord = params.cheb.order;
a = zeros(4*ord,1);

for i = 1:4
    a((i-1)*ord+1:i*ord) = U_1_cheb(:,i);
end

phi_cheb_int = intval(1)*y.a1';
Ad_N = chebDF_intval(phi_cheb_int,ord,params);
A_N = Ad_N^-1;
a_bar = intval(1)*a;


Y0 = computeY0(A_N,a_bar,phi_cheb_int,params,ord,intICvec,params.del);
Y0hat = computeY0hat(A_N,rho,params.Lbvp,params.del,params.nu,a_bar(1:600),phi_cheb);
fprintf('Y bounds computed, Y0 = %d, Y0hat = %d\n', mid(Y0), mid(Y0hat));

Z0 = computeZ0(A_N,Ad_N,ord,params.del);
Z1 = computeZ1(A_N,ord,phi_cheb_int,params.del,params);
Z2hat = computeZ2hat(A_N,params.Lbvp,params.del,params.nu,rho);
fprintf('Z bounds computed, Z0 = %d, Z1 = %d, Z2hat = %d\n', mid(Z0), mid(Z1), mid(Z2hat));

rs = 0:10^-12:10^-6;
radii_poly = Y0 + Y0hat - (1-Z0-Z1-Z2hat)*rs;

good_r = sup((Y0 + Y0hat)/(1-Z0-Z1-Z2hat))

% pulse3 soln validated by leaving chebstar2 without fft in y3tail of Y0
% function (takes a while)
end