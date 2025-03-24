function [M_a_out] = M_op(a_in,w_s,mu_s)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

a0_prod = 4*( 1+(-mu_s)^2)* K_op(w_s,mu_s,1) + 4*K_op(w_s,mu_s,3);
a1_prod = 2*( 1+(-mu_s)^2)* w_s - 4* mu_s* K_op(w_s,mu_s,1)+6*K_op(w_s,mu_s,2);
a2_prod = -2*mu_s*w_s +4*K_op(w_s,mu_s,1)  ;

K1_a_in = K_op(a_in,mu_s,1);
K2_a_in = K_op(a_in,mu_s,2);

a0 = conv2(a_in,a0_prod);
a1 = conv2(K1_a_in,a1_prod);
a2 = conv2(K2_a_in,a2_prod);

M_a_out = a0+a1+a2;

end