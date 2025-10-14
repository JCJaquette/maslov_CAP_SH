clear all 

params.mu = 0.05; 
params.nu = 1.6; 
params.scale = 1/70;
params.lambda = 0;

order = 10-1;  % Manifold

params.order = order; % TODO Why do we have two orders?
params.mfld.order = order; 

% Add This boolean for intervals
params.isIntval = 1;

if params.isIntval
    params.mu = intval(params.mu);
    params.nu = intval(params.nu);
end


% -----------------------------------------------------------------------
% Now we calculate the coefficients of the parameterization for the stable
% and unstable manifold up to a desired order. 
 
% TODO: Replace convolution with Taylor_Convolution.

%  TODO: Make get_mflds universal
mflds=get_mflds(params);

D = diag([mflds.values.s,mflds.values.u]);
V = [mflds.vectors.s,mflds.vectors.u];


uncoeff   = mflds.unstable.coeffs;
stabcoeff = mflds.stable.coeffs;


figure(1)
tiledlayout(2,2) 

% Plotting stuff
uncoeff_plot = uncoeff;
stabcoeff_plot = stabcoeff;
if params.isIntval
    uncoeff_plot =uncoeff_plot.mid;
    stabcoeff_plot =stabcoeff_plot.mid;
end


nexttile
plot_coeff(uncoeff_plot,order);
title('Coefficient norms for unstable parameterization')

nexttile
plot_coeff(stabcoeff_plot,order);
title('Coefficient norms for stable parameterization')

nexttile
plot_manifold(uncoeff_plot,order,'r');
title('Unstable Manifold')

nexttile
plot_manifold(stabcoeff_plot,order,'b');
title('Stable Manifold')


%--------------------------------------------------------------------------
% Now we apply Lemma 4.4 to validate the parameterization we computed 

% TODO: Fix Intvals
 unpoly=mfld_poly(params,uncoeff,V,D);
 mflds.unstable.r_min = unpoly;
unstable_bound =unpoly;

 % stpoly=mfld_poly(params,stabcoeff,V,D);
 % mflds.stable.r_min = stpoly;
 % stable_bound = stpoly;
%  
  disp('The error on the parameterization for the unstable manifold is: ')
  disp(unstable_bound)
  % disp('The error on the parameterization for the stable manifold is: ')
 % disp(stable_bound)

 % TODO: Store the error bound

% Compute L_minus
 Lminus = computeLminus(params,mflds) ;
params.Lminus=Lminus;

% [vectors, values]= getJacEigs(0, params);
sigma_0 = exp(-real(mflds.values.u(1)) * params.Lminus)