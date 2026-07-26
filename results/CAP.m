clear 

%% Manually define parameters

Case_number = 1; 
% Case 1 : params.mu = 0.1; 
% Case 2 : params.mu = 0.1; 
% Case 3 : params.mu = 0.2;  

BOOL_load_bndl = 1;
bndl_BOOL.save_data = 0;
BOOL_load_pulse = 1;
BOOL_save_pulse = 0;
BOOL_load_Euminus = 1;
BOOL_save_Euminus = 0;

% Parameters for the pulse validation
params.rho = .99; %This is rho from section 2 in paper 3
params.tol=4e-14;
params.bd_scale = .2; %This sets how close the pulse gets to the manifold when we cut it off
params.new = 1.01; %new=delta in the paper, nu in the code for pulse existence CAP


    %ODE Parameters
if Case_number == 1 || Case_number == 2 
    params.mu = 0.1; 
    params.scale = .3;
    params.order = 40; %manifold taylor coeffs
    params.pulse.order = 2^10; %cheb coeffs of pulse
    params.Eu.order = 2^9; %cheb coeffs of Eu-
else
    params.mu = 0.2; 
    params.scale = .3;
    params.order = 26; %manifold taylor coeffs
    params.pulse.order = 2^11; %cheb coeffs of pulse
    params.Eu.order = 2^10; %cheb coeffs of Eu-
end

if Case_number == 2 
    params.xi = pi;
else
    params.xi = 0;
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


%% Bundle and Manifold Computation

if BOOL_load_bndl 

    if Case_number == 1 || Case_number == 2 
        data_str = "data_bndl_nu_1p6_mu_0p1";
    else
        data_str = "data_bndl_nu_1p6_mu_0p2";
    end

    load(data_str)
    disp(['Loaded bundles and manifolds for case ',int2str(Case_number)])

    if Case_number == 2
        params.xi = pi;
    end
else
    
    % Computation 
    bndl_BOOL.plot = 0;
    bndl_BOOL.save_image = 0;
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

    disp('Computing pulse')

    % Get seed for Newton
    seed = get_newton_seed(params,mflds);     
    
    % Refine with Newton
    pulse4D = refine_cheb_orbit(seed,mflds,params);
    
    disp('Validating Pulse')

    % Validate
    [verif, pulse4D.r, vali_data] = verify_homoclinic_orbit(params,mflds,pulse4D,params.new);

    if BOOL_save_pulse
    
        if Case_number == 1
            data_str = "data_pulse1_nu_1p6_mu_0p1";
        elseif Case_number == 2 
            data_str = "data_pulse2_nu_1p6_mu_0p1";
        else
            data_str = "data_pulse_nu_1p6_mu_0p2";
        end

        save(data_str,'pulse4D')
        disp(['Saved pulse ',int2str(Case_number)])

    end

end

%    '_PO',int2str(params.pulse.order)];
% 
% z = flattenstruct(vali_data, '');
% data_Table = struct2table(z,'AsArray',true);
% 
% save(file_str,'data_Table')


%% Eu- Computation

if BOOL_load_Euminus

    if Case_number == 1
        data_str = "data_Eu1_nu_1p6_mu_0p1";
    elseif Case_number == 2 
        data_str = "data_Eu2_nu_1p6_mu_0p1";
    else
        data_str = "data_Eu_nu_1p6_mu_0p2";
    end

    load(data_str)

else

    disp('Computing Eu-')
    if Case_number == 3
        params.Eu.order = 2^10;   
    end
    params.del = 1.01;%two dels?
    %Get U_{\varphi'} and U_1
    Eu = chebInt(params,mflds,pulse4D);  

    disp('Getting CAP for Eu-')

    Eu.r = EuminusCAP(params,mflds,pulse4D,Eu);

    if BOOL_save_Euminus

        if Case_number == 1
            data_str = "data_Eu1_nu_1p6_mu_0p1";
        elseif Case_number == 2 
            data_str = "data_Eu2_nu_1p6_mu_0p1";
        else
            data_str = "data_Eu_nu_1p6_mu_0p2";
        end
    
        save(data_str,'params','Eu')
        disp(['Saved Eu- for case ',int2str(Case_number)])

    end

end

%% L+ Computation

disp('Computing L+')

mflds.Lplus = computeLplus(params,bndl,mflds,pulse4D,Eu.U_1_cheb);

%% Counting Zeros/Conjugate Points

disp('Counting zeros of determinant')

zerocount = 0;
zerofinder_tol = 1e-5;

BOOLzf_plot = 1;

disp('Finding zeros on [-L_conj,-L_bvp]')
zerocount = zerocount + countBeforeBVP(params,mflds,pulse4D,zerofinder_tol,BOOLzf_plot);

disp('Finding zeros on [-L_bvp,L_bvp]')
zerocount = zerocount + countBVP(pulse4D,Eu,zerofinder_tol,BOOLzf_plot);

disp('Finding zeros on [L_bvp,L_conj]')
zerocount = zerocount + countAfterBVP(params,bndl,mflds,pulse4D,Eu,zerofinder_tol,BOOLzf_plot);

