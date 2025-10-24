% @ N=32, with intlab 
% Takes ~13 min to get manifold 
% 6 min per manifold validation
% 20 sec to get bundles 
% 10 min to validate bundle
 
%% Initialize Parameters
clear

% ODE Parameters
params.mu = 0.2; 
params.nu = 1.6;
% Computational Parameters
params.scale = 1/200;
order = 2; 
% Interval Arithmetic
params.isIntval =1;
if params.isIntval 
    params.mu = intval('0.2'); 
    params.nu = intval('1.6');
end



% Computation 
BOOL_plot = 1;
BOOL_save_image = 0;
BOOL_save_data =0;
BOOL_Lminus = 0;

% Potential parameter for finding 
params.lambda = 0; 

% Setting several things in memory
if params.isIntval
    zero=intval(0);
else
    zero=0;
end
params.order = order; 
params.mfld.order = order; 


%% Get Manifolds
tic
%  TODO: Make get_mflds universal
mflds=get_mflds(params);
time_get_mflds = toc 


tic
BOOL_stable = 1;
disp('Computing Radii Poly Bounds')
[mflds,r_min_s,data_mfld_poly_s ]=mfld_poly(params, mflds,BOOL_stable );
r_min_s

 time_mfld_poly = toc

if isnan(r_min_s )
    return
end
% return

%% Compute L_minus
 if BOOL_Lminus 
    BOOL_stable = 0;
    [mflds,r_min_u,data_mfld_poly_u]=mfld_poly(params, mflds,BOOL_stable );

    Lminus = computeLminus(params,mflds) ;
    params.Lminus=Lminus;
    sigma_0 = exp(-real(mflds.values.u(1)) * params.Lminus)
 end
% return

%% Get Bundles
tic
% Computes Manifold coeff, and bundle Coeff.
[bndl] = getAllBundleCoefficients(params,mflds);

time_get_bndl = toc 

tic
disp('Computing Radii Poly Bounds')
[ r_min, data_bndl_poly] = bundle_rad_poly(params,mflds,bndl);
time_bndl_poly = toc
r_min


%% Plot

% Plot manifold and bundles
if BOOL_plot 
    %% Plot Bundles
    figure
    plots=plot_bndl(params,mflds,bndl,'b');
    if BOOL_save_image
        obj= gca;
        exportgraphics(obj,'manifold_bndl.png',Resolution=500)
    end


    %% Plot Coefficients
    figure
    % TODO 
    plot_coeff_sum(mflds.stable.coeffs,params,'o');
    hold on 
    bundle_permute = permute(bndl.coeffs,[3,4,1,2]);
    bundle_stab = reshape(bundle_permute(:,:,:,1),[order+1,order+1,4]);
    plot_coeff_sum((bundle_stab) ,params,'^');
    % plot_coeff_sum(imag(bundle_stab) ,params)
    bundle_unstab = reshape(bundle_permute(:,:,:,3),[order+1,order+1,4]);
    plot_coeff_sum((bundle_unstab ) ,params,'square');
    % plot_coeff_sum(imag(bundle_unstab ) ,params)

    if params.isIntval
        plot_coeff_sum((intval(mflds.stable.coeffs.rad) ) ,params,'.');
        plot_coeff_sum((intval(bundle_stab.rad) ) ,params,'*');
        plot_coeff_sum((intval(bundle_unstab.rad) ) ,params,'x');
    end

    legend('stable manifold','stable bundle','unstable bundle','error manifold','error s bundle','error u bundle')
    xlabel('$n$','Interpreter','latex')
    xlim([0,order])
    ylim([-25,2])
    
    if BOOL_save_image
        obj= gca;
        exportgraphics(obj,'Coeff_size.png',Resolution=500)
    end
end
if BOOL_save_data 
    save('nu_1p6_mu_0p2_V2')
end
 
return
