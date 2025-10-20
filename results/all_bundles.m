% Takes ~11 min to run with intlab
 
%% Initialize Parameters
clear
x = zeros(1, 4); 
params.lambda = 0; 
params.mu = 0.2; 
params.nu = 1.6; 
params.scale = 1/2;
params.isIntval =0;

BOOL_plot = 1;
BOOL_save = 0;
BOOL_Lminus = 1;

if params.isIntval
    zero=intval(0);
else
    zero=0;
end

order = 32; 
params.order = order; 
params.mfld.order = order; 


%% Get Manifolds
tic
%  TODO: Make get_mflds universal
mflds=get_mflds(params);
time_get_mflds = toc 
tic
BOOL_stable = 1;
[mflds,r_min_s]=mfld_poly(params, mflds,BOOL_stable );
r_min_s
 time_mfld_poly = toc


%% Compute L_minus
 if BOOL_Lminus 
    BOOL_stable = 0;
    [mflds,r_min_u]=mfld_poly(params, mflds,BOOL_stable );

    Lminus = computeLminus(params,mflds) ;
    params.Lminus=Lminus;
    sigma_0 = exp(-real(mflds.values.u(1)) * params.Lminus)
 end


%% Get Bundles
tic
% Computes Manifold coeff, and bundle Coeff.
[bndl] = getAllBundleCoefficients(params,mflds);

time_get_bndl = toc 

tic
disp('Computing Radii Poly Bounds')
[ r_min ] = bundle_rad_poly(params,mflds,bndl);
time_bndl_poly = toc
r_min


%% Plot

% Plot manifold and bundles
if BOOL_plot 
    %% Plot Bundles
    figure
    plots=plot_bndl(params,mflds,bndl,'b');
    if BOOL_save
        obj= gca;
        exportgraphics(obj,'manifold_bndl.png',Resolution=500)
    end


    %% Plot Coefficients
    figure
    plot_coeff_sum(mflds.stable.coeffs,params,'o')
    hold on 
    bundle_permute = permute(bndl.coeffs,[3,4,1,2]);
    bundle_stab = reshape(bundle_permute(:,:,:,1),[order+1,order+1,4]);
    plot_coeff_sum((bundle_stab) ,params,'x')
    % plot_coeff_sum(imag(bundle_stab) ,params)
    bundle_unstab = reshape(bundle_permute(:,:,:,3),[order+1,order+1,4]);
    plot_coeff_sum((bundle_unstab ) ,params,'square')
    % plot_coeff_sum(imag(bundle_unstab ) ,params)
    legend('stable manifold','stable bundle','unstable bundle')
    xlabel('$n$','Interpreter','latex')
    xlim([0,order])
    ylim([-16,2])
    
    if BOOL_save
        obj= gca;
        exportgraphics(obj,'Coeff_size.png',Resolution=500)
    end
end

 
return