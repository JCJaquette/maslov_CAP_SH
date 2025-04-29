function [ Y_0_out] = calc_Y_0(A_N,sol,mu_s,vs1_shift,vs2_shift,g_1,params)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

% Get the index;
N = sqrt(length(sol)+1)-1;
mu=params(1);

fout = F_op(sol,mu_s,vs1_shift,vs2_shift,g_1,params);

[N_out,M_out]=size(fout);
if N_out~=M_out
    disp("WARNING NON SQUARE DIMENSION!")
end

f_N = fout(1:N+1,1:N+1);

f_N_vec =  reshape(f_N,(N+1)^2,1);
f_N_vec(1)=[];

Af_N = A_N*f_N_vec;

f_infty=fout;
f_infty(1:N+1,1:N+1) = 0*fout(1:N+1,1:N+1);

tic
[ L_entry,~ ] = L_op(N_out-1,mu_s,mu);
toc 
% Need to prevent dividing by zero
L_entry(1:3,1:3)=L_entry(1:3,1:3)+pi^3;

Af_infty = f_infty./L_entry;
Af_infty(1:N+1,1:N+1) = 0*Af_infty(1:N+1,1:N+1);

Y_0_out = sum(abs(Af_N),"all") + sum(abs(Af_infty ),"all");




end