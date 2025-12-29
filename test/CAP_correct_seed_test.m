clear %update to github
%varbs for this is sec1

%To load the variables I was working with there are p1sec1 and p3sec1, 
% which gives the variables that would come from just running the first 
% section of the code with pulse1 and pulse3 respectively.

%ODE Parameters
params.mu = 0.05; 
params.nu = 1.6;

% For pulse 3(mu,nu = .2,1.6): scale = .3, order = 15

% Computational Parameters
params.scale = .09;
params.order = 15; 
params.mfld.order = params.order;

% Interval Arithmetic
params.isIntval =0;
if params.isIntval 
    params.mu = intval(params.mu); 
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
[mflds,bndl,Lminus] = get_all_bundles(params,bndl_BOOL); 

% NOTE: Removed "mflds_r" from get_all_bundles output. This data is stored
% in "mflds.stable.r_min" or "mflds.unstable.r_min" 
% Also removed "bndl_r"; this is stored in the bndl object

% TODO: This should not be recast as a double at the top level. 
% If you need to recast it as a double, do it inside the necessary function. 
[mflds,intradii] = struct_intvaltodouble(mflds);


% params.stable.error = mflds_r + max(intradii.stable.coeffs(1,2,:)); this
% isn't final, maybe not exactly right?
% params.unstable.error = mflds_r + max(intradii.unstable.coeffs(1,2,:));

return

%%
% sec2

% Parameters for the pulse validation
params.rho = .99;
params.cheb.order=2^9;
params.tol=4e-14;
params.Lbvp = 0;
params.bd_scale = .2;
params.new = 1.01;
params.xi = 0;

% TODO: This should not be recast as a double at the top level. 
% If you need to recast it as a double, do it inside the necessary function. 
params = struct_intvaltodouble(params);

% We have three pulses(xi is the branch):
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

%%

%plot what goes into newton vs what comes out for this pulse 

figure
hold on
x = linspace(-1,1,201);
seedFN = chebSum(seed.a1',-1:.01:1);
plot(x,seedFN)
newtFN = chebSum(y.a1',-1:.01:1);
plot(x,newtFN)

%%

%do the same for the old psoln
%not exactly the same as in the old files since mflds could be different
%things don't seem to work as well here, maybe go back to old files on github?

% TODO: This should not be recast as a double at the top level. 
% If you need to recast it as a double, do it inside the necessary function. 
params = struct_intvaltodouble(params);
load('psoln1.mat')
soln = psoln1;

[seeed,params.Lbvp] = old_get_newton_seed(soln,params,mflds);

y1 = refine_cheb_orbit(seeed,mflds,params);
figure
hold on
x = linspace(-1,1,201);
seeedFN = chebSum(seeed.a1',-1:.01:1);
plot(x,seeedFN)
newtFN1 = chebSum(y1.a1',-1:.01:1);
plot(x,newtFN1)


%%

% Validate
verify_homoclinic_orbit(params,mflds,y,params.new);

%%

load('psoln1.mat')
plot(psoln1(:,1),psoln1(:,5))
hold on
plot(psoln(:,1),psoln(:,5))
