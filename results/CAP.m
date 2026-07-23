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
return
disp('Counting zeros of determinant')

zerocount = 0;
zerofinder_tol = 10e-3;
BOOLzf_plot = 1;

disp('Finding zeros on [-L_conj,-L_bvp]')
zerocount = zerocount + countBeforeBVP(params,mflds,pulse4D,zerofinder_tol,BOOLzf_plot);

disp('Finding zeros on [-L_bvp,L_bvp]')
zerocount = zerocount + countBVP(pulse4D,Eu,zerofinder_tol,BOOLzf_plot);

disp('Finding zeros on [L_bvp,L_conj]')
zerocount = zerocount + countAfterBVP(params,bndl,mflds,pulse4D,Eu,zerofinder_tol,BOOLzf_plot);

%%

if 1
    
    S = [1, 0, 0, 0; 
         0, 0, 1, 0;
         0, 2, 0, 1;
         0, 1, 0, 0];

    [tbeta,tgamma] = getLinSolnCoords(params,bndl,Eu.U_1_cheb,pulse4D);
    teta = getLinSolnCoords(params,bndl,Eu.U_vpp_cheb,pulse4D);
    tildes = [tbeta,tgamma,teta];

    ax_pts=intval([]);bx_pts=ax_pts;V1s_pts=ax_pts;V2s_pts=ax_pts;
    domain = linspace(0,mflds.Lplus,1000);

    for x = domain

        [ax_pts(:,end+1),bx_pts(:,end+1),V1s_pts(:,end+1),V2s_pts(:,end+1)] = get_a_b(params,bndl,mflds,pulse4D,tildes,x,S);

    end

    dadx_ode = ODEtest(params,pulse4D,mflds,ax_pts,domain);
    dbdx_ode = ODEtest(params,pulse4D,mflds,bx_pts,domain);
    dV1dx_ode = ODEtest(params,pulse4D,mflds,V1s_pts,domain);
    dV2dx_ode = ODEtest(params,pulse4D,mflds,V2s_pts,domain);

    dx = domain(2) - domain(1);

    for i = 2:length(domain(2:end-1))

        dadx_fd(:,i) = (ax_pts(:,i+1) - ax_pts(:,i-1))/(2*dx);
        dbdx_fd(:,i) = (bx_pts(:,i+1) - bx_pts(:,i-1))/(2*dx);
        dV1dx_fd(:,i) = (V1s_pts(:,i+1) - V1s_pts(:,i-1))/(2*dx);
        dV2dx_fd(:,i) = (V2s_pts(:,i+1) - V2s_pts(:,i-1))/(2*dx);

    end


figure
tiledlayout(2,2)

nexttile
plot(domain,mid(dadx_ode(1,:)))
hold on
plot(domain(2:end-1),mid(dadx_fd(1,:)))
title('da1/dx')

nexttile
plot(domain,mid(dbdx_ode(1,:)))
hold on
plot(domain(2:end-1),mid(dbdx_fd(1,:)))
title('db1/dx')

nexttile
plot(domain,mid(dV1dx_ode(1,:)))
hold on
plot(domain(2:end-1),mid(dV1dx_fd(1,:)))
title('dV1s/dx')

nexttile
plot(domain,mid(dV2dx_ode(1,:)))
hold on
plot(domain(2:end-1),mid(dV2dx_fd(1,:)))
title('dV2s/dx')  

sgtitle('domain(2:end-1)')

figure
tiledlayout(2,2)

nexttile
plot(domain,mid(dadx_ode(1,:)))
hold on
plot(domain(1:end-2),mid(dadx_fd(1,:)))
title('da1/dx')

nexttile
plot(domain,mid(dbdx_ode(1,:)))
hold on
plot(domain(1:end-2),mid(dbdx_fd(1,:)))
title('db1/dx')

nexttile
plot(domain,mid(dV1dx_ode(1,:)))
hold on
plot(domain(1:end-2),mid(dV1dx_fd(1,:)))
title('dV1s/dx')

nexttile
plot(domain,mid(dV2dx_ode(1,:)))
hold on
plot(domain(1:end-2),mid(dV2dx_fd(1,:)))
title('dV2s/dx') 

sgtitle('domain(1:end-2)')

end

function [ax,bx,V1s,V2s] = get_a_b(params,bndl,mflds,pulse4D,tildes,x,S)

    tbeta = tildes(:,1);
    tgamma = tildes(:,2);
    teta = tildes(:,3);

    sig = get_sig_afterBVP(params,pulse4D,x);
    sig1 = sig(1); sig2 = sig(2);
    get_Vs(params,bndl,mflds,sig(1),sig(2),x,S);

    tV1s = [exp(mflds.values.s(1) *x);0;0;0];
    tV2s = [0;exp(mflds.values.s(2) *x);0;0];
    
    a_res = [bndl.normalForm(1,3,3,1), bndl.normalForm(1,4,2,2);
             bndl.normalForm(2,3,2,2), bndl.normalForm(2,4,1,3)];
    
    tV1u = [a_res(1,1)*x * sig1^2 * exp(mflds.values.s(1) *x);
            a_res(2,1)*x * sig1*sig2 * exp(mflds.values.s(2) *x);
            exp(mflds.values.u(1) *x);
            0];
    tV2u = [a_res(1,2)*x * sig1*sig2 * exp(mflds.values.s(1) *x);
            a_res(2,2)*x * sig2^2 * exp(mflds.values.s(2) *x);
            0;
            exp(mflds.values.u(2) *x)];
    
    W_sig = bndl_one_point(real(sig1),imag(sig1),bndl,params);
    W_sig = W_sig + infsup(-bndl.r_min.sup,bndl.r_min.sup);

    V1s = S * (W_sig*tV1s);
    V2s = S * (W_sig*tV2s);
    V1u = S * (W_sig*tV1u);
    V2u = S * (W_sig*tV2u);

    ax = teta(1)*V1s + teta(2)*V2s;
    
    bx = tbeta(1)*V1s + tbeta(2)*V2s + tgamma(1)*V1u + tgamma(2)*V2u;
    
    ax = real(ax);
    bx = real(bx);

end

function dfdx = ODEtest(params,pulse4D,mflds,fx,xrange)

nu = params.nu; mu = params.mu;

    for i = 1:length(xrange)

        sig = get_sig_afterBVP(params,pulse4D,xrange(i));

        Px = mfld_one_point(real(sig(1)),imag(sig(1)),mflds.stable.coeffs,params);
        w1 = Px(1);

        dfdx(:,i) = [0,0,0,1;
                     0,0,1,-2;
                     -(1+mu) + 2*nu*w1 - 3*w1^2, 0,0,0;
                     0,1,0,0]*fx(:,i);

    end

end
