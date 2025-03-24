function [resonant_bundle_strc] = Compute_Resonant_Bundles(w_s,mu_s)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

resonant_bundle_strc =0;


[N,M] = size(w_s);
N=max([N,M]);
N=N-1;

N_inv = 2*N;

[ K_out ] = K_op(w_s,mu_s,1);


% Define operator L & M & K

% Define w^I which is the approximate convolution inverse of w_s
[w_I] = inverse_convolution(w_s,N_inv);



% A = Cauchy_Mat_rep(w_I,N)
% A = Cauchy_Mat_rep(w_s,2*N)
% A = Cauchy_Mat_rep_trunc(w_s,N_inv)
% Ai=inv(A);
% w_ii = Ai(:,1);
% w_i_vec = reshape(w_ii,N_inv+1,N_inv+1)
% 
% 23;
% return
% x=A*w_s_vec
% y=conv2(w_s,w_I)
% xx=reshape(x,N_inv+N+1,N_inv+N+1)
% return
% w_I_mat = A(:,1);
% 
% 
% w_I_mat =reshape(w_I_mat,N+1,N+1)
% w_I-w_I_mat
% 
% return
[LinvM_mat] = LinvM(w_I,w_s,mu_s);

B_dagger_N_N = eye((N+1)^2)+LinvM_mat;
 
B=inv(B_dagger_N_N);
% Define Bdagger ( as a matrix )

% Define B_{N,N} 
% Define B_{infty,N} 
% Define A 

% define T- newton operator
% Define approximate solution. 

% keyboard


end
