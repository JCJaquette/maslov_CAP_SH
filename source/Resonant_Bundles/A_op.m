function [A_out] = A_op(w_I,w_s,B_N_N,mu_s,a_in,N)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

% A= B * \hat{L}^{-1}

wI_a = conv2(a_in,w_I);

La = K_op(wI_a,mu_s,-3);

La_pi_N = La(1:N+1,1:N+1);
La_pi_infty = La;
La_pi_infty(1:N+1,1:N+1) =0*La_pi_infty(1:N+1,1:N+1);

La_size = size(La_pi_infty);
N_La = La_size(1);

La_N_vec = reshape(La_pi_N,(N+1)^2,1);
% La_N_vec(1)=[]; 
% B_La_N= 0*(1:(N+1)^2);
% B_La_N(2:end) = B_N_N*La_N_vec;
B_La_N = B_N_N*La_N_vec;
B_La_N = reshape(B_La_N,N+1,N+1);

% compute -Bdagger_{infty,N}

MBLa = M_op(B_La_N,w_s,mu_s);

LiMBLa = conv2(MBLa,w_I);

LiMBLa = K_op(LiMBLa,mu_s,-3);

LiMBLa(1:N+1,1:N+1)=0*LiMBLa(1:N+1,1:N+1);

LiMBLa_size = size(LiMBLa);
N_LiMBLa = LiMBLa_size(1);

N_infty = max([N_LiMBLa,N_La]);
A_out=zeros(N_infty,N_infty);

% Add the infty part
A_out(1:N_LiMBLa,1:N_LiMBLa)= -LiMBLa;
A_out(1:N_La,1:N_La)=A_out(1:N_La,1:N_La)+La_pi_infty;
% A_out=0*A_out;
% add the finite part
A_out(1:N+1,1:N+1)=B_La_N;

end