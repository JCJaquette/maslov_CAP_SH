function [resonant_bundle_strc] = Compute_Resonant_Bundles(w_s,mu_s)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

resonant_bundle_strc =0;

w_s=0*w_s;
w_s(1,1)=1;
w_s(1,2)=1e-5;
w_s(2,1)=1e-5;

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
% Includes 0,0 entry
[LinvM_mat] = LinvM(w_I,w_s,mu_s);

B_dagger_N_N = eye((N+1)^2)+LinvM_mat;
B_dagger_N_N(1,1)=0;

B_dagger_N_N_restrict = B_dagger_N_N(2:end,2:end);
B_N_N_restrict = inv(B_dagger_N_N_restrict);

[V0,D0]=eig(B_dagger_N_N_restrict);
D_diag0= diag(D0);
[~,ord0]=sort(abs(D_diag0));
D_diag0=D_diag0(ord0);
V_sort0 = V0(:,ord0);

vec_0 = V_sort0(:,1);
vec_0=[0;vec_0];
vec_0 = reshape(vec_0,1+N,1+N);

[V,D]=eig(B_dagger_N_N);
ddd=diag(D);
[~,ord]=sort(abs(ddd));
ddd=ddd(ord);

% V_sort = V(:,ord);
Inver_E = 1./D_diag0;
Inver_E(1)=0;
Approx_inv = V_sort0*diag(Inver_E)*inv(V_sort0);

xxxx= B_dagger_N_N_restrict*Approx_inv;

% d0 =diag(D0);
%  hold on 
% scatter(real(D_diag0),imag(D_diag0))
% scatter(real(ddd),imag(ddd))

% keyboard

% TODO the matrix B_N_N is horrible ill conditioned. 

% invert thing except for the zero eval
[V,D]=eig(B_dagger_N_N);
D_diag= diag(D);
[~,ord]=sort(abs(D_diag));
V_sort = V(:,ord);
Inver_E = 1./D_diag(ord);
Inver_E(1)=0;
Inver_E(2)=0;

vec_0 = V_sort(:,1);
vec_0 = reshape(vec_0,1+N,1+N);

vec_1 = V_sort(:,2);
vec_1 = reshape(vec_1,1+N,1+N);

% Approx_inv = V_sort*diag(Inver_E)*inv(V_sort);

figure
scatter(real(D_diag),imag(D_diag))
figure

B_N_N_old=inv(B_dagger_N_N);
B_N_N=Approx_inv;

e_00=0*w_s;
e_00(1,1)=1;


f_0 = M_op(e_00,w_s,mu_s);
f_0=f_0(1:N+1,1:N+1);

N_ApproxMax = 5*N+1;
w_u_bar = zeros(N_ApproxMax ,N_ApproxMax );

% Af0_out = A_op(w_I,w_s,B_N_N_restrict,mu_s,-f_0,N);


% w_u_bar = Af0_out;
% w_u_bar(1,1)=1; 

vn2 = sum(abs(vec_0.*conj(vec_0)),'all');

disp('Improving')
for i = 1:20
        disp('Run next')
    [F_out] = F_op(w_u_bar,w_s,mu_s);


    sum(abs(F_out),'all')

    Y0 = A_op(w_I,w_s,B_N_N,mu_s,F_out,N);
    % Y0 = A_op(w_I,w_s,B_N_N,mu_s,F_out,N);
    sum(abs(Y0),'all')

    w_u_bar = w_u_bar-Y0(1:N_ApproxMax ,1:N_ApproxMax );
    
    % Project away from kernel
    % inner = sum(w_u_bar(1:1+N,1:1+N).*conj(vec_0),'all');

    % w_u_bar(1:1+N,1:1+N) = w_u_bar(1:1+N,1:1+N) - (inner/vn2)*vec_0;
    % sum(w_u_bar(1:1+N,1:1+N).*conj(vec_0),'all')

    surf(log(abs(F_out))/log(10))
    % surf(log(abs(Y0))/log(10))
%     xlim([0,30])
% ylim([0,30])
% zlim([-50,10])
    keyboard
end


% sum(abs(F_out),'all')

% approx_sol(1,1)=-1; 
% [F_out] = F_op(approx_sol,w_s,mu_s);

% sum(abs(F_out),'all')


zlim([-50,10])


keyboard
% Define Bdagger ( as a matrix )

% Define B

% Define B_{N,N} 
% Define B_{infty,N} 
% Define A 

% define T- newton operator
% Define approximate solution. 

% keyboard


end
