function [ K_out ] = K_op(a_in,mu_s,power)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here
% Computes K^s a_in
% mu_s - eigenvalue

[N,M]=size(a_in);
if N~=M
    disp("WARNING NON SQUARE DIMENSION!")
end
N=N-1;

diff_vec = 0:N;

diff_mat = zeros(N+1,N+1);
for i=1:N+1
    diff_mat(1:end,i) = diff_vec';
end


left = diff_mat*mu_s;
K_entry = left+left';
% raise to power
K_entry =(K_entry.^power);
K_entry(1,1)=0;

K_out = K_entry.*a_in;

end