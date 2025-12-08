clear %update to github

%ODE Parameters
params.mu = 0.2; 
params.nu = 1.6;

% For pulse 3(mu,nu = .2,1.6): scale = .3, order = 15

% Computational Parameters
params.scale = .3;
params.order = 15; 
params.mfld.order = params.order;

% Interval Arithmetic
params.isIntval =1;
if params.isIntval 
    params.mu = intval('0.2'); 
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

% Computation 
bndl_BOOL.plot = 0;
bndl_BOOL.save_image = 0;
bndl_BOOL.save_data = 0;
bndl_BOOL.Lminus = 1; 
bndl_BOOL.stable = 1;

% Get the bundles and manifolds
[mflds,mflds_r,bndl,bndl_r,Lminus] = all_bundles(params,bndl_BOOL); 

[mflds,intradii] = struct_intvaltodouble(mflds);
params.stable.error = mflds_r + max(intradii.stable.coeffs(1,2,:));
params.unstable.error = mflds_r + max(intradii.unstable.coeffs(1,2,:));

%%

% Parameters for the pulse validation
params.rho = .99;
params.cheb.order=2^10;
params.tol=4e-14;
params.Lbvp = 0;
params.bd_scale = .1;
params.new = 1.01;
params.xi = 0;

params = struct_intvaltodouble(params);

% We have three pulses:
% mu=.05,nu=1.6, xi=0 
% mu=.05,nu=1.6, xi=π
% mu=.2,nu=1.6, xi=0

% Get seed for Newton
[seed,params.Lbvp] = get_newton_seed(params,mflds);


% Refine with Newton
y = refine_cheb_orbit(seed,mflds,params);

% yo1 = chebcoeff_to_function(new_y.a1);
% yo2 = chebcoeff_to_function(new_y.a2);
% yo3 = chebcoeff_to_function(new_y.a3);
% yo4 = chebcoeff_to_function(new_y.a4);
% plot_manifold(mflds.unstable.coeffs, params.mfld.order, 'red')
% plot_manifold(mflds.stable.coeffs, params.mfld.order, 'blue')
% plot3(yo1 ,yo2 ,yo4 ,'LineWidth',1)

hold on
seedFN = chebSum(seed.a1',-1:.01:1);
plot(seedFN)
newtFN = chebSum(y.a1',-1:.01:1);
plot(newtFN)

% Validate
verify_homoclinic_orbit(params,mflds,y,params.new);


%%


chebInt %This gets U_{\varphi'} and U_1
clearvars -except mflds bndl new_y phiPrime_cheb h_cheb 

% load('saved_stuff/sig0_1.mat'); %Pick sig0 depending on the n
% load('saved_stuff/sig0_2.mat');
load('saved_stuff/sig0_3.mat');

computeLplus(params,bndl,mflds,U_1,sig0);
