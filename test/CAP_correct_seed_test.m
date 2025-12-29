clear 

CASE_number = 3; 
% Case 1 : params.mu = 0.05; 
% Case 2 : params.mu = 0.05; 
% Case 3 : params.mu = 0.20; 

BOOL_load_bndl =1 ;

% Debugging Ideas:
%  @@ Adjust the getting of the pulse, so that it uses fsolve first like Hannah did
%  #### It seems that this matches up. .... So continue with the debugging!
%  
%  @@ Add some quantifiable tests that output how close we are to being
%  successful!
%  @@ Next, test to see if the boundary conditions, and L_bvp are different

% Load bundles or recompute
if BOOL_load_bndl 
    if CASE_number == 1 || CASE_number == 2 
        data_str = "data_bndl_nu_1p6_mu_0p05";
    else
        data_str = "data_bndl_nu_1p6_mu_0p2";
    end
    load(data_str)
else
    % Computational Parameters
    params.scale = .09;
    params.order = 15; 
    params.mfld.order = params.order;
    params.isIntval =0;
    
    % Computation 
    bndl_BOOL.plot = 0;
    bndl_BOOL.save_image = 0;
    bndl_BOOL.save_data = 1;
    bndl_BOOL.Lminus = 1; 
    bndl_BOOL.stable = 1;

    
    %ODE Parameters
    if CASE_number == 1 || CASE_number == 2 
        params.mu = 0.05; 
    else
        params.mu = 0.2; 
    end
    params.mu = 0.2; 
    params.nu = 1.6;
    
    % For pulse 3(mu,nu = .2,1.6): scale = .3, order = 15
    
    
    % Interval Arithmetic    
    if params.isIntval 
        params.mu = intval(params.mu);  % This line doesn't faithfull cast as interval enclosure of '.2'
        params.nu = intval('1.6');
    end
    
    % Potential parameter for finding 
    params.lambda = 0; 
    
    % % Setting several things in memory
    % if params.isIntval
    %     zero=intval(0);
    % else
    %     zero=0;
    % end
    

    
    % Get the bundles and manifolds
    [mflds,bndl] = get_all_bundles(params,bndl_BOOL); 
    
    % NOTE: Removed "mflds_r" from get_all_bundles output. This data is stored
    % in "mflds.stable.r_min" or "mflds.unstable.r_min" 
    % Also removed "bndl_r"; this is now stored in the bndl object
    % Also removed "Lminus"; this is now stored in the mflds object
    
    % TODO: This should not be recast as a double at the top level. 
    % If you need to recast it as a double, do it inside the necessary function. 
    [mflds,intradii] = struct_intvaltodouble(mflds);
    
    
    
    % params.stable.error = mflds_r + max(intradii.stable.coeffs(1,2,:)); this
    % isn't final, maybe not exactly right?
    % params.unstable.error = mflds_r + max(intradii.unstable.coeffs(1,2,:));

    % NOTE: This ↑↑↑ error calculation should be done where it is used, not here
    % Also, if you are not immediately sure, then this is something
    % nontrivial, that merits working out on paper / latex
end



%%
% sec2

% Parameters for the pulse validation
params.rho = .99;
params.cheb.order=2^9;
params.tol=4e-14;
params.Lbvp = 0;
params.bd_scale = .2;
params.new = 1.01;
% We have three pulses(xi is the branch):
if CASE_number == 2 
    params.xi = pi;
else
    params.xi = 0;
end
% mu=.05,nu=1.6, xi=0 
% mu=.05,nu=1.6, xi=π
% mu=.2,nu=1.6, xi=0


params = struct_intvaltodouble(params);
% TODO: This ↑↑↑ should not be recast as a double at the top level!!!! 
% If you need to recast it as a double, do it inside the necessary function. 



% Get seed for Newton
[seed,params.Lbvp] = get_newton_seed(params,mflds);
 title('New seed')

% Refine with Newton
y = refine_cheb_orbit(seed,mflds,params);

% yo1 = chebcoeff_to_function(new_y.a1);
% yo2 = chebcoeff_to_function(new_y.a2);
% yo3 = chebcoeff_to_function(new_y.a3);
% yo4 = chebcoeff_to_function(new_y.a4);
% plot_manifold(mflds.unstable.coeffs, params.mfld.order, 'red')
% plot_manifold(mflds.stable.coeffs, params.mfld.order, 'blue')
% plot3(yo1 ,yo2 ,yo4 ,'LineWidth',1)

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
 
% verify_homoclinic_orbit(params,mflds,y,params.new);



%%


figure
plot(soln(:,1),soln(:,5))
title('Initial, full solution')
% hold on
% plot(psoln(:,1),psoln(:,5))
