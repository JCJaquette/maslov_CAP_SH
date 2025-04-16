function [resonant_bundle_strc] = Compute_Resonant_Bundles_2(mu_s,vs1,g_1,params)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

resonant_bundle_strc =0;
[N,M] = size(vs1);
N=max([N,M]);
N=N-1;

vs2=vs1';

vs1_shift = -4*(K_op(vs1,mu_s,3)+K_op(vs1,mu_s,1));
vs2_shift = -4*(K_op(vs2,mu_s,3)+K_op(vs2,mu_s,1));
% Shift Down
shift1_down = diag(1+0*(1:N),-1);
vs1_shift=shift1_down*vs1_shift;
vs2_shift=shift1_down*vs2_shift;
% K_op(vs1,mu_s,3)

[A_dag] = A_dagger(mu_s,g_1,vs1_shift,vs2_shift,params);
null(A_dag)
A_dag_short=A_dag(2:end,2:end);
disp('null')
null(A_dag_short)

g_1_vec = reshape(g_1,(N+1)^2,1);

g_1_vec(1)=[];

sol = A_dag_short\g_1_vec;

sol_mat = [0;sol];
sol_mat =reshape(sol_mat,N+1,N+1);


surf(log(abs(sol_mat))/log(10))

keyboard
return 

 

 
[LinvM_mat] = LinvM(w_I,w_s,mu_s);


end
