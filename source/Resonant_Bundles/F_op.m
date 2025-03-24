function [F_out] = F_op(a_in,w_s,mu_s)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

a_in(1,1)=1;

K3_a_in = K_op(a_in,mu_s,3);
La = conv2(w_s,K3_a_in);

M_a_out = M_op(a_in,w_s,mu_s);

F_out = La + M_a_out ;

end