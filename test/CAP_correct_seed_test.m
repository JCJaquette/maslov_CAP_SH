clear 
figure 
CASE_number = 3; 
% Case 1 : params.mu = 0.05; 
% Case 2 : params.mu = 0.05; 
% Case 3 : params.mu = 0.20; 

BOOL_load_bndl =0;

BOOL_validate_pulse = 1; 

    % Computational Parameters
    params.isIntval =1;
    
    % Computation 
    bndl_BOOL.plot = 0;
    bndl_BOOL.save_image = 0;
    bndl_BOOL.save_data = 1;
    bndl_BOOL.Lminus = 1; 
    bndl_BOOL.stable = 1;

    % It is still somewhat strange how the Chebyshev seed and Converge are
    % so different. There is still probably something weird with how
    % computational parameters are being chosen.
    % 
    % 
%
%  @@@@@ Mike, also, please add more comments to the files! What I'd recommend doing, at least, 
%       is while you are debugging something, add comments to the thing you're debugging. 
%       This helps with double checking that what you did in the code is correct 

 

% Load bundles or recompute
if BOOL_load_bndl 
    if CASE_number == 1 || CASE_number == 2 
        % data_str = "data_bndl_nu_1p6_mu_0p05";
        data_str = "data_bndl_nu_1p6_mu_0p1";
    else
        data_str = "data_bndl_nu_1p6_mu_0p2";
    end
    load(data_str)
else
    
    %ODE Parameters
    if CASE_number == 1 || CASE_number == 2 
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
    
end

return

%%
% sec2

% Parameters for the pulse validation
params.rho = .99;%This is delta_s in paper 3
params.cheb.order=2^10;
params.tol=4e-14;
params.Lbvp = 0;
params.bd_scale = .2;%This sets how close the pulse gets to the manifold when we cut it off
params.new = 1.01;%new=delta in the paper, nu in the code for pulse existence CAP
% We have three pulses(xi is the branch):
if CASE_number == 2 
    params.xi = pi;
else
    params.xi = 0;
end
% mu=.05,nu=1.6, xi=0 
% mu=.05,nu=1.6, xi=π
% mu=.2,nu=1.6, xi=0


BOOL_load_oldpsoln = 0;

if BOOL_load_oldpsoln
    params.loadpsoln = Case_number;
else
    params.loadpsoln = 0;
end

% Get seed for Newton
[seed,params.Lbvp] = get_newton_seed(params,mflds);
 title('New seed')

% Refine with Newton
y = refine_cheb_orbit(seed,mflds,params);



%%

%plot what goes into newton vs what comes out for this pulse 

figure
hold on
x = linspace(-1,1,201);
seedFN = chebSum(seed.a1',-1:.01:1);
plot(x,seedFN)
newtFN = chebSum(y.a1',-1:.01:1);
plot(x,newtFN)
hold off
grid on
legend('New seed','Chebyshev Converged')
title('New seed, and convergant')

% Validate stuff with the new pulse, loaded from data
if BOOL_validate_pulse 
    verify_homoclinic_orbit(params,mflds,y,params.new);
end

return
% input('Continue?')
% return
%%

%do the same for the old psoln
%not exactly the same as in the old files since mflds could be different
%things don't seem to work as well here, maybe go back to old files on github?


params = struct_intvaltodouble(params);
% TODO: This ↑↑↑ should NOT be recast as a double at the top level. 
% If you need to recast it as a double, do it inside the necessary function. 

if CASE_number == 1 
    load('psoln1.mat')
    soln = psoln1;
elseif CASE_number == 2 
    load('psoln2.mat')
    soln = psoln2;
else
    load('psoln3.mat')
    soln = psoln3;
end


figure
[seeed,params.Lbvp] = old_get_newton_seed(soln,params,mflds);
title('Old seed')

y1 = refine_cheb_orbit(seeed,mflds,params);

figure
hold on
x = linspace(-1,1,201);
seeedFN = chebSum(seeed.a1',-1:.01:1);
plot(x,seeedFN)
newtFN1 = chebSum(y1.a1',-1:.01:1);
plot(x,newtFN1)
hold off
grid on 
legend('Old seed','Chebyshev Converged')
title('Old seed and convergant')

%%

figure 
title('Two Seeds')
hold on

plot(x,seedFN)
plot(x,seeedFN)
legend('new','old')
grid on 
hold off

figure 
title('Two Convergants')

hold on
plot(x,newtFN)
plot(x,newtFN1)
legend('new','old')
grid on 
hold off

%%

% Validate
 
% Validate stuff with the old pulse, loaded from data
if BOOL_validate_pulse 
    verify_homoclinic_orbit(params,mflds,y1,params.new);
end



%%


figure
plot(soln(:,1),soln(:,5))
title('Initial, full solution')
% hold on
% plot(psoln(:,1),psoln(:,5))
