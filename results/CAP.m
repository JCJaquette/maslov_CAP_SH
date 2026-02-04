clear 

%% Manually define parameters

Case_number = 3; 
% Case 1 : params.mu = 0.1; 
% Case 2 : params.mu = 0.1; 
% Case 3 : params.mu = 0.2;  

BOOL_load_bndl = 1;
BOOL_load_pulse = 1;
BOOL_load_Euminus = 0;

% Parameters for the pulse validation
params.rho = .99; %This is delta_s in paper 3
params.pulse.order = 2^10; %cheb coeffs of pulse
params.Eu.order = 2^9; %cheb coeffs of Eu-
params.tol=4e-14;
params.bd_scale = .2;%This sets how close the pulse gets to the manifold when we cut it off
params.new = 1.01;%new=delta in the paper, nu in the code for pulse existence CAP


    %ODE Parameters
if Case_number == 1 || Case_number == 2 
    params.mu = 0.1; 
    params.scale = .3;
    params.order = 40;
else
    params.mu = 0.2; 
    params.scale = .3;
    params.order = 26;  
end
params.nu = 1.6;
params.mfld.order = params.order;
  
% Interval Arithmetic  
params.isIntval = 1;
if params.isIntval 
    params.mu = intval(num2str(params.mu));  
    params.nu = intval('1.6');
end

% Potential parameter for finding 
params.lambda = 0; 

%% Automatically define parameters

if Case_number == 2 
    params.xi = pi;
else
    params.xi = 0;
end

%% Bundle and Manifold Computation

if BOOL_load_bndl 

    if Case_number == 1 || Case_number == 2 
        data_str = "data_bndl_nu_1p6_mu_0p1";
    else
        data_str = "data_bndl_nu_1p6_mu_0p2";
    end

    load(data_str)
    disp(['Loaded bundles and manifolds for case ',int2str(Case_number)])

else
    
    % Computation 
    bndl_BOOL.plot = 0;
    bndl_BOOL.save_image = 0;
    bndl_BOOL.save_data = 1;
    bndl_BOOL.Lminus = 1; 
    bndl_BOOL.stable = 1;
    
    % For pulse 3(mu,nu = .2,1.6): scale = .3, order = 15
    
    
    % Get the bundles and manifolds
    [mflds,bndl] = get_all_bundles(params,bndl_BOOL); 

end


%% Pulse Computation


if BOOL_load_pulse

    if Case_number == 1
        data_str = "data_pulse1_nu_1p6_mu_0p1";
    elseif Case_number == 2 
        data_str = "data_pulse2_nu_1p6_mu_0p1";
    else
        data_str = "data_pulse_nu_1p6_mu_0p2";
    end

    load(data_str)
    disp(['Loaded pulse ',int2str(Case_number)])

else    
    
    % We have three pulses(xi is the branch):
    % mu=.1,nu=1.6, xi=0 
    % mu=.1,nu=1.6, xi=π TODO:this w/ mu=.1
    % mu=.2,nu=1.6, xi=0

    disp('Computing pulse')

    % Get seed for Newton
    [seed,params.Lbvp] = get_newton_seed(params,mflds);     
    
    % Refine with Newton
    pulse4D = refine_cheb_orbit(seed,mflds,params);
    
    disp('Validating Pulse')

    % Validate
    [verif, pulse4D.r] = verify_homoclinic_orbit(params,mflds,pulse4D,params.new);

end

%% Eu- Computation

if BOOL_load_Euminus

    if Case_number == 1
        data_str = "data_pulse1_nu_1p6_mu_0p1";
    elseif Case_number == 2 
        data_str = "data_pulse2_nu_1p6_mu_0p1";
    else
        data_str = "data_pulse_nu_1p6_mu_0p2";
    end

    load(data_str)

else

    disp('Computing Eu-')
    
    params.del = 1.01;%two dels?
    
    %Get U_{\varphi'} and U_1
    [U_vp, U_1, params.nonzero] = chebInt(params,mflds,pulse4D); 
    %We may work with lower order here (but still keep some tail of 0s for proof)
    phi_cheb = pulse4D.a1(1:params.nonzero)';  
    
    params.rho = get_rho(params,mflds,pulse4D,phi_cheb);
    
    r = EuminusCAP(params,phi_cheb,U_1);

end

%% L+ Computation
disp('Computing L+')

computeLplus(params,bndl,mflds,U_1,sig0);
