clear all 

%-------------------------------------------------------------------------
% first we calculate the size of the eigenvector/eigenvalue enclosures
% !!! TODO: Make Rigorous; Maybe use formula for e-vec?
% 
point=[0,0,0,0];


params.mu = 0.05; 
params.nu = 1.6; 
params.scale = 1/7;
params.lambda = 0;

order = 10-1;  % Manifold

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
 
%  TODO: Make get_mflds universal
mflds=get_mflds(params);

D = diag([mflds.values.s,mflds.values.u]);
V = [mflds.vectors.s,mflds.vectors.u];


uncoeff   = mflds.unstable.coeffs;
stabcoeff = mflds.stable.coeffs;


k=6;
mat=zeros(k^2,4)*stabcoeff(1,1,1);
K = zeros(k^2, 2);
for i = 1:k 
    for j = 1:k
        row=(i-1)*k + j;
        mat(row,:)=stabcoeff(i,j,:);
        K(row,:)=[i-1,j-1];
    end  
end

figure
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
 % return
 unstable_bound = min(unpoly(find(unpoly > 0)));

 stpoly=mfld_poly(params,stabcoeff,V,D);
 stable_bound = min(stpoly(find(stpoly > 0)));
%  
  disp('The error on the parameterization for the unstable manifold is: ')
  disp(unstable_bound)
  disp('The error on the parameterization for the stable manifold is: ')
 disp(stable_bound)

 % TODO: Store the error bound

% Compute L_minus
 Lminus = computeLminus(params,mflds) ;
params.Lminus=Lminus;

% [vectors, values]= getJacEigs(0, params);
sigma_0 = exp(-real(mflds.values.u(1)) * params.Lminus)