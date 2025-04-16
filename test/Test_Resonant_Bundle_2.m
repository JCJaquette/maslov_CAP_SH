clear

cd ..
cd("results")
load ('Manifold_DATA.mat')
% load ('Manifold_DATA_50.mat')
cd ..
cd("test")
 
[N,n,j] = size(stabcoeff);
N=N-1;

mu_s = stabeigs(1);
mu_s2 = stabeigs(1);

% Differentiate with respect to the first index
diff = 0:N;
diff_mat = diag(diff);

% Only interested in the first coefficient
non_res_coeff = stabcoeff(:,:,1);

mu=params(1);
nu=params(2);

g_11 = -3*conv2(non_res_coeff,non_res_coeff);
g_11(1:N+1,1:N+1)= 2*nu*non_res_coeff + g_11(1:N+1,1:N+1);
g_1= g_11(1:N+1,1:N+1);
% TODO Imncrease size.

vs1 = diff_mat*non_res_coeff ;
% vs2 = non_res_coeff*diff_mat;
vs2=vs1';

Compute_Resonant_Bundles_2(mu_s,vs1,g_1,params)

return 
% Shift up
non_res_coeff(1:end-1,:)=non_res_coeff(2:end,:);
non_res_coeff(end,:)=0*non_res_coeff(end,:);

w_s = non_res_coeff(:,:,1);
w_s=w_s/w_s(1,1);
% Trim
% www=w_s;
% w_s=zeros(2*N+1,2*N+1);
% w_s(1:N+1,1:1+N)=www;
% w_s(:,end)=[];
% w_s(end,:)=[];



mu_s = stabeigs(1);

Compute_Resonant_Bundles(w_s,mu_s)
 

