function [LinvM_mat] = LinvM(w_I,w_s,mu_s)
%UNTITLED8 Summary of this function goes here
%   Detailed explanation goes here

[n,m]=size(w_s)
beta = real(mu_s);
gamma = imag(mu_s);
N=n-1;
a0_prod = 4*( 1+(-mu_s)^2)* K_op(w_s,mu_s,1) + K_op(w_s,mu_s,3);
a1_prod = ( 1+(-mu_s)^2)* w_s - 4* mu_s* K_op(w_s,mu_s,1)+6*K_op(w_s,mu_s,2);
a2_prod = -2*mu_s*w_s +4*K_op(w_s,mu_s,1)  ;

a0 = conv2(w_I,a0_prod);
a1 = conv2(w_I,a1_prod);
a2 = conv2(w_I,a2_prod);

[N1,N2]=size(a0);

G_op = zeros(N1,N1);

for n1=0:N1-1
    for n2=0:N1-1
        if n1+n2==0
            val =1;
        elseif beta^2/gamma^2 >= abs(n1-n2)/(n1+n2)
            val =1;
        else
            num=(n1+n2)^2*beta^2+(n1-n2)^2*gamma^2;
            num = sqrt((beta^2+gamma^2)*num);
            den = 2 *beta*gamma*max([n1,n2]);
            val = num/den;
        end
        G_op(n1+1,n2+1)=val;
    end
end

bound = (beta^2+gamma^2)/(2*beta*gamma);

(1/(beta*N))^3*sum(abs(a0),'all')
(1/(beta*N))^2*sum(abs(a1),'all')*bound
(1/(beta*N))^1*sum(abs(a2),'all')*bound^2

(1/(beta*N))^3*sum(abs(a0),'all')
(1/(beta*N))^2*sum(G_op.*abs(a1),'all')
(1/(beta*N))^1*sum((G_op.^2).*abs(a2),'all')

% Truncate ;; This needs to be fixed 
a0=a0(1:N+1,1:N+1);
a1=a1(1:N+1,1:N+1);
a2=a2(1:N+1,1:N+1);


% K_entry 

diff_vec = 0:N;

diff_mat = zeros(N+1,N+1);
for i=1:N+1
    diff_mat(1:end,i) = diff_vec';
end
left = diff_mat*mu_s;
K_entry = left+left';

% % Kinv_entry 
% % 
% % diff_vec = 0:2*N+1;
% % diff_mat = zeros(2*N+2,2*N+2);
% % for i=1:2*N+2
% %     diff_mat(1:end,i) = diff_vec';
% % end
% % left = diff_mat*mu_s;
% % K_entry_inv = 1./(left+left'); 
K_entry_inv =1./K_entry
K_entry_inv(1,1)=0;

K_mat=reshape(K_entry,(N+1)^2,1);
K_inv3_mat=reshape(K_entry_inv.^3,(N+1)^2,1);

K_mat=diag(K_mat);
K_inv3_mat=diag(K_inv3_mat);
%%%

LinvM_mat_0 = K_inv3_mat*Cauchy_Mat_rep_trunc(a0,N );
LinvM_mat_1 = K_inv3_mat*Cauchy_Mat_rep_trunc(a1 ,N)*K_mat;
LinvM_mat_2 = K_inv3_mat*Cauchy_Mat_rep_trunc(a2 ,N)*K_mat*K_mat;

LinvM_mat=LinvM_mat_0+LinvM_mat_1+LinvM_mat_2;

end