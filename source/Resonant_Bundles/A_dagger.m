function [A_dag] = A_dagger(mu_s,g_1,vs1_shift,vs2_shift,params)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
A_dag=0;

[n,m]=size(vs1_shift)
beta = real(mu_s);
gamma = imag(mu_s);
N=n-1;

mu=params(1);

g_mat = Cauchy_Mat_rep_trunc(g_1,N );
 
% Leading Operator
L_entry = zeros(N+1,N+1);
for m=0:N
    for n=0:N
        xi = -(1+((m-1)*mu_s+n*mu_s')^2)^2-mu;
        if (m==2&&n==0)||(m==1&&n==1)
            xi=0; % Going to be zero anywhay
        end
        L_entry(m+1,n+1)=xi;
    end
end

L_mat=reshape(L_entry,(N+1)^2,1); 


L_mat=diag(L_mat); 


%  Note! Tricky how this is organized. 
% vs1 (2,0) % 3
% vs2 (1,1) % (N+1)+2
vs1_vec = reshape(vs1_shift,(N+1)^2,1);
vs2_vec = reshape(vs2_shift,(N+1)^2,1);

ResSol = zeros((N+1)^2,(N+1)^2);
ResSol(:,3) = vs1_vec;
ResSol(:,(N+1)+2) = vs2_vec;

g_mat(:,3)=0* g_mat(:,3);
g_mat(:,(N+1)+2)=0* g_mat(:,(N+1)+2);
A_dag = L_mat+g_mat+ResSol;
rank(A_dag);

return

 
end