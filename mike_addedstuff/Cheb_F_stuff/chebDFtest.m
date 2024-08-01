clear %copy pasted from pulse validation to get cheb coeffs
% Set initial seed for Newton's method 
[x,params, mflds] = init_sol_val_tests();

params.normalizeBasis = 0;

disp('Choice of Parameters')
disp(params)
disp(['Chebyshev order: ', num2str(params.cheb.order)])
disp(['Manifold order: ', num2str(params.mfld.order)])

disp(['Initial choice for phi1, phi2, psi: ', num2str(x.phi1), ', ', num2str(x.phi2), ', ', num2str(x.psi)]);


u_ic = [-0.190343484893489 ; -0.128234936941383 ; 0.193875175006471 ; 0.168311653800752];
s_ic = [-0.190343484893489 ; 0.128234936941383 ; 0.193875175006471 ; -0.168311653800752];

L1=params.L;


[V, D]=asym_un_vecs(params);
v1=real(V(:,1));

%%% Option to use NF solution
sol=BK_nf_4dim(params,0,L1, 0);
time_vec=sol(:,1);
sol=sol(:,2:5);


cheb_coeff = get_cheb_coeffs(sol, params);
x.a1=cheb_coeff(:,1)';
x.a2=cheb_coeff(:,2)';
x.a3=cheb_coeff(:,3)';
x.a4=cheb_coeff(:,4)';

nozeros=abs(nonzeros(x.a1));
max_ord=my_order(max(nozeros));
min_ord=my_order(min(nozeros));

disp(['Decay in Chebyshev coefficients for first component: ', num2str(max_ord-min_ord)])

G=FHomoclinic(x, mflds, params);

% Perform Newton's method 

new_x = refine_cheb_orbit(x, mflds, params);
% new_x.a1 gives b

phi_cheb = new_x.a1;
%%
% making sure my DF function works by comparing with numerical derivative
% at random point in random direction

h1 = rand(1,1800);
h2 = rand(1,1800);

eps = 1e-8;


fperth = (chebF(h1 +eps*h2,2*[1,1,1,1],phi_cheb,450,params) - chebF(h1,2*[1,1,1,1],phi_cheb,450,params))/eps;

dfh = chebDF(phi_cheb,450,params)*h2';

diff = fperth - dfh;

norm(diff)







