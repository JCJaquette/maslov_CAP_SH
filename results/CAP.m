clear %update to github
%varbs for this is sec1

%To load the variables I was working with there are p1sec1 and p3sec1, 
% which gives the variables that would come from just running the first 
% section of the code with pulse1 and pulse3 respectively.

Case_number = 3; 
% Case 1 : params.mu = 0.1; 
% Case 2 : params.mu = 0.1; 
% Case 3 : params.mu = 0.2;  

BOOL_load_bndl =1;

BOOL_load_pulse = 1;

if BOOL_load_bndl 

    if Case_number == 1 || Case_number == 2 
        % data_str = "data_bndl_nu_1p6_mu_0p05";
        data_str = "data_bndl_nu_1p6_mu_0p1";
    else
        data_str = "data_bndl_nu_1p6_mu_0p2";
    end

    load(data_str)

else
    % Computational Parameters
    params.isIntval = 1;
    
    % Computation 
    bndl_BOOL.plot = 0;
    bndl_BOOL.save_image = 0;
    bndl_BOOL.save_data = 1;
    bndl_BOOL.Lminus = 1; 
    bndl_BOOL.stable = 1;

    
    %ODE Parameters
    if Case_number == 1 || Case_number == 2 
        params.mu = 0.1; 
        params.scale = .3;
        params.order = 20;
    else
        params.mu = 0.2; 
        params.scale = .3;
        params.order = 26;  
        params.cheb.order=2^10;
    end
    params.nu = 1.6;
    params.mfld.order = params.order;
    
    % For pulse 3(mu,nu = .2,1.6): scale = .3, order = 15
    
    
    % Interval Arithmetic    
    if params.isIntval 
        params.mu = intval(num2str(params.mu));  
        params.nu = intval('1.6');
    end
    
    % Potential parameter for finding 
    params.lambda = 0; 
    
    
    % Get the bundles and manifolds
    [mflds,bndl] = get_all_bundles(params,bndl_BOOL); 
    
    % NOTE: Removed "mflds_r" from get_all_bundles output. This data is stored
    % in "mflds.stable.r_min" or "mflds.unstable.r_min" 
    % Also removed "bndl_r"; this is now stored in the bndl object
    % Also removed "Lminus"; this is now stored in the mflds object

end

return
%%
% sec2


if BOOL_load_pulse

    if Case_number == 1
        % data_str = "data_bndl_nu_1p6_mu_0p05";
        data_str = "data_pulse1_nu_1p6_mu_0p1";
    elseif Case_number == 2 
        data_str = "data_pulse2_nu_1p6_mu_0p1";
    else
        data_str = "data_pulse_nu_1p6_mu_0p2";
    end

    load(data_str)

else

    % Parameters for the pulse validation
    params.rho = .99;%This is delta_s in paper 3
    params.cheb.order=2^10;
    params.tol=4e-14;
    params.Lbvp = 0;
    params.bd_scale = .2;%This sets how close the pulse gets to the manifold when we cut it off
    params.new = 1.01;%new=delta in the paper, nu in the code for pulse existence CAP
    params.xi = 0;
    
    
    % We have three pulses(xi is the branch):
    % mu=.1,nu=1.6, xi=0 
    % mu=.05,nu=1.6, xi=π TODO:this w/ mu=.1?
    % mu=.2,nu=1.6, xi=0
    
    % Get seed for Newton
    [seed,params.Lbvp] = get_newton_seed(params,mflds);
     
    
    % Refine with Newton
    pulse4D = refine_cheb_orbit(seed,mflds,params);
    
    % yo1 = chebcoeff_to_function(new_y.a1);
    % yo2 = chebcoeff_to_function(new_y.a2);
    % yo3 = chebcoeff_to_function(new_y.a3);
    % yo4 = chebcoeff_to_function(new_y.a4);
    % plot_manifold(mflds.unstable.coeffs, params.mfld.order, 'red')
    % plot_manifold(mflds.stable.coeffs, params.mfld.order, 'blue')
    % plot3(yo1 ,yo2 ,yo4 ,'LineWidth',1)
    
    figure
    hold on
    x = linspace(-params.Lbvp,params.Lbvp,201);
    seedFN = chebSum(seed.a1',-1:.01:1);
    plot(x,seedFN)
    newtFN = chebSum(pulse4D.a1',-1:.01:1);
    plot(x,newtFN)
    legend('into newton','out of newton')
    xlabel('t')
    ylabel('$\varphi$(t)',Interpreter='latex')
    
    % Validate
    verify_homoclinic_orbit(params,mflds,pulse4D,params.new);

end

%%

params.cheb.order = 600;
params.del = 1.01;
% TODO: get this stuff manually vv
maxphi = max(abs(pulse4D.phi1),abs(pulse4D.phi2)); 
mani_error = 3.2e-11; 
ICerror = 2*pi/log(1/maxphi) * mani_error;
params.rho = max(2.4e-8,ICerror);
%                               ^^
manifold_u.coeffs = mflds.unstable.coeffs;


[U_vp, U_1] = chebInt(params,mflds,pulse4D); %This gets U_{\varphi'} and U_1

EuminusCAP(params,pulse4D,U_1)


%%

% load('saved_stuff/sig0_1.mat'); %Pick sig0 depending on the n
% load('saved_stuff/sig0_2.mat');
load('saved_stuff/sig0_3.mat');

computeLplus(params,bndl,mflds,U_1,sig0);
