clear all 

params.mu = 0.2; 
params.nu = 1.6; 
params.scale = 1/2;
params.lambda = 0;

order = 30;  % Manifold

params.order = order; % TODO Why do we have two orders?
params.mfld.order = order; 

% Add This boolean for intervals
params.isIntval = 0;

if params.isIntval
    params.mu = intval(params.mu);
    params.nu = intval(params.nu);
end


% -----------------------------------------------------------------------
% Now we calculate the coefficients of the parameterization for the stable
% and unstable manifold up to a desired order. 
 
% TODO: Replace convolution with Taylor_Convolution. ??

%  TODO: Make get_mflds universal
mflds=get_mflds(params);

 


figure(1)
tiledlayout(2,2) 

% Plotting stuff
uncoeff_plot = mflds.unstable.coeffs;
stabcoeff_plot = mflds.stable.coeffs;
if params.isIntval
    uncoeff_plot =uncoeff_plot.mid;
    stabcoeff_plot =stabcoeff_plot.mid;
end


nexttile
plot_coeff(uncoeff_plot,order);
title('Coefficient norms for unstable parameterization')

nexttile
% plot_coeff(stabcoeff_plot,order);
plot_coeff_sum(mflds.stable.coeffs,params,'o')
title('Coefficient norms for stable parameterization')

nexttile
plot_manifold(uncoeff_plot,order,'r');
title('Unstable Manifold')

nexttile
plot_manifold(stabcoeff_plot,order,'b');
title('Stable Manifold')

figure(20)
plot_manifold(uncoeff_plot,order,'r');
hold on
plot_manifold(stabcoeff_plot,order,'b');
hold off
%--------------------------------------------------------------------------
% Now we apply Lemma 4.4 to validate the parameterization we computed 



% TODO: Fix Intvals
BOOL_stable = 0;
 [mflds,r_min_u] =mfld_poly(params,mflds,BOOL_stable);
 % mflds.unstable.r_min = unstable_bound;
 
BOOL_stable =1;
 [mflds,r_min_s]=mfld_poly(params,mflds,BOOL_stable );
 % mflds.stable.r_min = stable_bound;
 
  disp('The error on the parameterization for the unstable manifold is: ')
  disp(r_min_u)
  disp('The error on the parameterization for the stable manifold is: ')
 disp(r_min_s)

 % TODO: Store the error bound

% Compute L_minus
 Lminus = computeLminus(params,mflds) ;
params.Lminus=Lminus;

% [vectors, values]= getJacEigs(0, params);
sigma_0 = exp(-real(mflds.values.u(1)) * params.Lminus)