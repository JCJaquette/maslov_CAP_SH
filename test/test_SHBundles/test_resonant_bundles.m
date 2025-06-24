clear
x = zeros(1, 4); 
params.lambda = 0; 
params.mu = 0.05; 
params.nu = 1.6; 
params.scale = 1/5;

order = 40-1; 
params.order = order; 
params.mfld.order = order; 
[eigenvectors, eigenvalues] = getJacEigs_toMerge(0, params); 

params.eigenvalues.s = eigenvalues.s; 
params.eigenvectors.s = eigenvectors.s; 
params.eigenvalues.u = eigenvalues.u; 
params.eigenvectors.u = eigenvectors.u; 


% Computes Manifold coeff, and bundle Coeff.
[All_Bundle_coeffs, normalForm_coeff] = getAllBundleCoefficients(params);



% Only nonzero components of normalForm_coeff are in 
% normalForm_coeff(i1,i2,i3,i4) for i1=1,2,  i2 =3,4, and i3+i4 =4.
% This is the tensor A, from the resonant bundle section.

% Plots norm of coefficietns
res_bundle_norm = sqrt(sum(abs(All_Bundle_coeffs).^2,1));
res_bundle_norm = permute ( res_bundle_norm, [3,4,2,1]);
res_bundle_norm_single_list = zeros(order+1,4);
for i=1:order+1
    local_norm =zeros(i,4);
    for j=1:i
        local_norm(j,:)=res_bundle_norm(j,i-j+1,:,1);
    end
    res_bundle_norm_single_list(i,:)=sum(local_norm,1);
end
plot(0:order,log(res_bundle_norm_single_list)/log(10),'o')
title('Log_{10} norm of the bundle coefficients')
legend('stab 1','stab 2','unst 1','unst 2')

% Picks out all of the stable / unstable bundles, and rearranges the data
% structure. 

stab_1 = reshape(All_Bundle_coeffs(:,1,:,:),[4 ,order + 1, order + 1]);
stab_2 = reshape(All_Bundle_coeffs(:,2,:,:),[4 ,order + 1, order + 1]);
unst_1 = reshape(All_Bundle_coeffs(:,3,:,:),[4 ,order + 1, order + 1]);
unst_2 = reshape(All_Bundle_coeffs(:,4,:,:),[4 ,order + 1, order + 1]);


stab_1  = permute ( stab_1, [2 3 1]);
stab_2  = permute ( stab_2, [2 3 1]);
unst_1 = permute ( unst_1, [2 3 1]);
unst_2 = permute ( unst_2, [2 3 1]);


return










