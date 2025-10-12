clear

tic 

x = zeros(1, 4); 
params.lambda = 0; 
params.mu = 0.05; 
params.nu = 1.6; 
params.scale = 1/7;
params.isIntval = 0;

if params.isIntval
    zero=intval(0);
else
    zero=0;
end

order = 20-1; 
params.order = order; 
params.mfld.order = order; 



% TODO: What is the difference between getJacEigs and getJacEigs_toMerge??
% params.eigenvalues = eigenvalues; 
% params.eigenvectors = eigenvectors; 

%  TODO: Make get_mflds universal
mflds=get_mflds(params);


% Computes Manifold coeff, and bundle Coeff.
[All_Bundle_coeffs, normalForm_coeff,mflds] = getAllBundleCoefficients(params,mflds);

time1 = toc 
tic


% 1. Get Manifold Coefficients
manifold_coeff = mflds.stable.coeffs;
manifold_coeff_norm=zeros(order+1,1)*zero;
% 2. Get Bundle Coefficients

disp('Computing Radii Poly Bounds')
[ r_min ] = bundle_rad_poly(params,All_Bundle_coeffs, normalForm_coeff,mflds);
time1 
toc
r_min

try
    [y,f]=audioread('JobDone.mp3');
    sound(y,1.05*f)
end
return
%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot



% Only nonzero components of normalForm_coeff are in 
% normalForm_coeff(i1,i2,i3,i4) for i1=1,2,  i2 =3,4, and i3+i4 =4.
% This is the tensor A, from the resonant bundle section.

% Plots norm of coefficients
res_bundle_norm = sqrt(sum(abs(All_Bundle_coeffs).^2,1));
res_bundle_norm = permute ( res_bundle_norm, [3,4,2,1]);
res_bundle_norm_single_list = zeros(order+1,4)*zero;
for alpha=1:order+1
    % For each order
    local_norm =zeros(alpha,4)*zero;
    local_norm_manifold = zeros(alpha,1)*zero; 
    for j=1:alpha
        % for each combo of the order
        local_norm(j,:)=res_bundle_norm(j,alpha-j+1,:,1);
        local_norm_manifold(j)=sqrt(sum(abs(manifold_coeff(j,alpha-j+1,:)).^2,'all'));
    end
    res_bundle_norm_single_list(alpha,:)=sum(local_norm,1);
    manifold_coeff_norm(alpha)=sum(local_norm_manifold,1);
end

% Y_0b = sum(res_bundle_norm_single_list);
% need Epsilon infty & K bound

if params.isIntval
    manifold_coeff_norm=manifold_coeff_norm.sup;
    res_bundle_norm_single_list=res_bundle_norm_single_list.sup;
end

plot(0:order,log(manifold_coeff_norm)/log(10),'o')
hold on 
plot(0:order,log(res_bundle_norm_single_list)/log(10),'o')
hold off 
title('Log_{10} norm of the bundle coefficients')
legend('manifold','stab 1','stab 2','unst 1','unst 2')

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










