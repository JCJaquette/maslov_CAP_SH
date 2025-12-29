function [mflds,bndl,Lminus] = get_all_bundles(params,BOOL)
%% Get Manifolds

tic
%  TODO: Make get_mflds universal
mflds=get_mflds(params);
time_get_mflds = toc 


tic
disp('Computing Radii Poly Bounds')
[mflds,r_min_s,data_mfld_poly_s ]=mfld_poly(params, mflds,BOOL.stable );

time_mfld_poly = toc

if isnan(r_min_s )
    return
end
% return

%% Compute L_minus
 if BOOL.Lminus 
    BOOL.stable = 0;
    [mflds,r_min_u,data_mfld_poly_u]=mfld_poly(params, mflds,BOOL.stable );

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

bndl.r_min = r_min;


%% Plot

% Plot manifold and bundles
if BOOL.plot 
    %% Plot Bundles
    figure
    plots=plot_bndl(params,mflds,bndl,'b');
    if BOOL.save_image
        obj= gca;
        exportgraphics(obj,'manifold_bndl.png',Resolution=500)
    end


    %% Plot Coefficients
    figure
    % TODO 
    plot_coeff_sum(mflds.stable.coeffs,params,'o');
    hold on 
    bundle_permute = permute(bndl.coeffs,[3,4,1,2]);
    bundle_stab = reshape(bundle_permute(:,:,:,1),[params.order+1,params.order+1,4]);
    plot_coeff_sum((bundle_stab) ,params,'^');
    % plot_coeff_sum(imag(bundle_stab) ,params)
    bundle_unstab = reshape(bundle_permute(:,:,:,3),[params.order+1,params.order+1,4]);
    plot_coeff_sum((bundle_unstab ) ,params,'square');
    % plot_coeff_sum(imag(bundle_unstab ) ,params)

    if params.isIntval
        plot_coeff_sum((intval(mflds.stable.coeffs.rad) ) ,params,'.');
        plot_coeff_sum((intval(bundle_stab.rad) ) ,params,'*');
        plot_coeff_sum((intval(bundle_unstab.rad) ) ,params,'x');
    end

    legend('stable manifold','stable bundle','unstable bundle','error manifold','error s bundle','error u bundle')
    xlabel('$n$','Interpreter','latex')
    xlim([0,params.order])
    ylim([-25,2])
    
    if BOOL.save_image
        obj= gca;
        exportgraphics(obj,'Coeff_size.png',Resolution=500)
    end
end
if BOOL.save_data 
    save('nu_1p6_mu_0p2_V2')
end
 

end