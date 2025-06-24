clear
x = zeros(1, 4); 
params.lambda = 0; 
params.mu = 0.05; 
params.nu = 1.6; 
params.scale = 1/5;

order = 32-1; 
params.order = order; 
params.mfld.order = order; 
[eigenvectors, eigenvalues] = getJacEigs_toMerge(0, params); 

params.eigenvalues.s = eigenvalues.s; 
params.eigenvectors.s = eigenvectors.s; 
params.eigenvalues.u = eigenvalues.u; 
params.eigenvectors.u = eigenvectors.u; 

% st_coeffs = getStBundleCoefficients(params);
unst_coeffs = getUnstBundleCoefficients(params);

% stab_1 = reshape(unst_coeffs(:,1,:,:),[4 ,order + 1, order + 1]);
% stab_2 = reshape(unst_coeffs(:,2,:,:),[4 ,order + 1, order + 1]);
% unst_1 = reshape(unst_coeffs(:,3,:,:),[4 ,order + 1, order + 1]);
% unst_2 = reshape(unst_coeffs(:,4,:,:),[4 ,order + 1, order + 1]);
% 
% stab_1  = permute ( stab_1, [2 3 1]);
% stab_2  = permute ( stab_2, [2 3 1]);
% 
% unst_1 = permute ( unst_1, [2 3 1]);
% unst_2 = permute ( unst_2, [2 3 1]);
% 
% stab_A = (stab_1 +stab_2)/2;
% stab_B = (stab_1 -stab_2)/(2*1i);
% 
% unst_A = (unst_1 +unst_2)/2;
% unst_B = (unst_1 -unst_2)/(2*1i);
% 
% [~,plotpoints_As ]=plot_manifold(stab_A,order,'b');
% [~,plotpoints_Bs ]=plot_manifold(stab_B,order,'b');
% [~,plotpoints_A ]=plot_manifold(unst_A,order,'b');
% [~,plotpoints_B ]=plot_manifold(unst_B,order,'b');
% 
% 
% 
% p=30;
% X=zeros(30,30);
% Y=zeros(30,30);
% Z=zeros(30,30);
% r=linspace(0,1,p);
% theta=linspace(0,2*pi, p);
% for j=1:p
%     for k=1:p
%         ps1s2 = zeros(4,1);
%         s1=r(j)*cos(theta(k));
%         s2=r(j)*sin(theta(k));
%         X(j,k)=s1;
%         Y(j,k)=s2;
%     end
% end
% 
% figure
% 
% quiver3(X,Y,Z,plotpoints_A(:,:,1),plotpoints_A(:,:,2),plotpoints_A(:,:,3))
% hold on 
% quiver3(X,Y,Z,plotpoints_B(:,:,1),plotpoints_B(:,:,2),plotpoints_B(:,:,3))
% 
% 
% figure
% 
% frame_mat = zeros(p,p,4,4);
% frame_matI = zeros(p,p,4,4);
% frame_mat_norm = zeros(p,p);
% frame_matI_norm = zeros(p,p);
% frame_mat(:,:,:,1)= plotpoints_As;
% frame_mat(:,:,:,2)= plotpoints_Bs;
% frame_mat(:,:,:,3)= plotpoints_A;
% frame_mat(:,:,:,4)= plotpoints_B;
% 
% s_dot_s = zeros(p,p);
% u_dot_u = zeros(p,p);
% 
% for i = 1:p
%     for j=1:p
% 
% 
%         local_frame = reshape(frame_mat(i,j,:,:),[4,4]);
%         local_frame_I = inv(local_frame );
%         frame_matI(i,j,:,:)=local_frame_I ;
% 
%         frame_mat_norm(i,j)= norm(local_frame);
%         frame_matI_norm(i,j)= norm(local_frame_I );
% 
%         s1=local_frame(:,1);
%         s2=local_frame(:,2);
%         u1=local_frame(:,3);
%         u2=local_frame(:,4);
% 
%         s_dot_s(i,j)=dot(s1,s2)/(norm(s1)*norm(s2));
%         u_dot_u(i,j)=dot(u1,u2)/(norm(u1)*norm(u2));
% 
%     end
% end
% 
% surf(X,Y,frame_matI_norm)
% 
% figure 
% surf(X,Y,s_dot_s)
% figure
% surf(X,Y,u_dot_u)


return
v_coeff_pm = reshape(st_coeffs(1,1,:,:), [order + 1, order + 1]); 
v_coeff_dp = bundle_coeff_from_st_mfld_deriv(params, 1);

growth_rate = eigenvalues.s(1);
x_range = .0001:.1:15; 

[x_range, u_dp] = scalar_sol_to_4D_sol(params, v_coeff_dp, growth_rate, x_range);
[x_range, u_pm] = scalar_sol_to_4D_sol(params, v_coeff_pm, growth_rate, x_range); 

u_dp = real(u_dp); 
u_pm = real(u_pm);

dv_coeff_dp = diff_non_res_variational_sol_coeff(params, v_coeff_dp, growth_rate, 1);
dv_coeff_pm = diff_non_res_variational_sol_coeff(params, v_coeff_pm, growth_rate, 1);

%% TODO Think about real vs complex valued quantities 
[new_x_range, dot_u_dp_c] = scalar_sol_to_4D_sol(params, dv_coeff_dp, growth_rate, x_range); 
[new_x_range, dot_u_pm] = scalar_sol_to_4D_sol(params, dv_coeff_dp, growth_rate, x_range); 

dot_u_dp = real(dot_u_dp_c); 
u_dp = real(u_dp);

params.mfld.st_coeffs = calc_proj_coeff(eigenvalues.s, eigenvectors.s, params); 
f_u_dp = SH_variational_eq(params, x_range, u_dp);
f_u_pm = SH_variational_eq(params, x_range, u_pm); 

should_be_0_dp = dot_u_dp - f_u_dp'; 
should_be_0_pm = dot_u_pm - f_u_pm'; 

disp("ERROR")
disp(vecnorm(should_be_0_pm))
disp(vecnorm(should_be_0_dp))


figure 
hold on 
plot(x_range, real(u_dp(:,1))*params.scale)
plot(x_range, real(u_pm(:,1)))












