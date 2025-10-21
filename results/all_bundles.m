% @ N=32, with intlab 
% Takes ~13 min to get manifold 
% 6 min per manifold validation
% 20 sec to get bundles 
% 10 min to validate bundle
 
%% Initialize Parameters
clear
x = zeros(1, 4); 
params.lambda = 0; 
params.mu = 0.2; 
params.nu = 1.6; 
params.scale = 1/2;
params.isIntval =1;

BOOL_plot = 1;
BOOL_save = 0;
BOOL_Lminus = 1;

if params.isIntval
    zero=intval(0);
else
    zero=0;
end

order = 30; 
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

% all_bundles
% time_get_mflds =
%      7.731630329000001e+02
% Computing Radii Poly Bounds
% Calculating Y0.
%   1.0e-016 *
%    0.7755481335____
% Calculating Z1.
%    0.3296918836028_
% Calculating Z2.
%     @(r)K_N*(6*linsum+2*params.nu+3*r)
%    0.7208189855424_
%    0.30130822692909
% r_min_s =
%      1.157014020859463e-16
% time_mfld_poly =
%      3.508071083000001e+02
% Calculating Y0.
%   1.0e-016 *
%    0.7755481335____
% Calculating Z1.
%    0.3296918836028_
% Calculating Z2.
%     @(r)K_N*(6*linsum+2*params.nu+3*r)
%    0.7208189855424_
%    0.30130822692909
% Found L_minus.
%   34.680107189078__
% intval sigma_0 = 
%   1.0e-003 *
%    0.5125868115903_
% Computing Manifold
% Computing Bundles
% time_get_bndl =
%   13.549899100000002
% Computing Radii Poly Bounds
% intval Ya_0 = 
%   1.0e-008 *
%    1.______________
% intval Yb_0 = 
%   1.0e-013 *
%    0.18240796______
% intval Yc_0 = 
%   1.0e-008 *
%    1.______________
% intval K_N = 
%    0.14305010031898
% intval Za = 
%    0.761361386189__
% intval Zb = 
%   1.0e-014 *
%    0.12478305441991
% intval Zc = 
%    0.1216362038152_
% intval Y_0 = 
%   1.0e-008 *
% [   0.00000182407958,   0.51792892341696] 
% intval Z = 
%    0.882997590004__
% intval r_min = 
%   1.0e-007 *
% [   0.00000155901026,   0.44266517539042] 
% time_bndl_poly =
%      5.762073050000001e+02
% intval r_min = 
%   1.0e-007 *
% [   0.00000155901026,   0.44266517539042] 
% Unable to perform assignment because value of type 'intval' is not convertible to 'double'.
% Error in plot_bndl (line 38)
%             plotpoints(j,k,:)=real(pt_local);
% Error in all_bundles (line 72)
%     plots=plot_bndl(params,mflds,bndl,'b');
% Caused by:
%     Error using double
%     Conversion to double from intval is not possible. 