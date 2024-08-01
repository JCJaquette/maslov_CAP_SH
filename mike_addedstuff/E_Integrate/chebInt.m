clear 
% copy pasted from pulse validation to get cheb coeffs of phi
% I changed the function init_sol_val_tests to get the evecs/evals out
% Set initial seed for Newton's method 
[x,params, mflds,evecs,evals] = init_sol_val_tests();

params.normalizeBasis = 0;

disp('Choice of Parameters')
disp(params)
disp(['Chebyshev order: ', num2str(params.cheb.order)])
disp(['Manifold order: ', num2str(params.mfld.order)])

disp(['Initial choice for phi1, phi2, psi: ', num2str(x.phi1), ', ', num2str(x.phi2), ', ', num2str(x.psi)]);


L1=params.L;


%%% Option to use NF solution
sol=BK_nf_4dim(params,0,L1, 0);
time_vec=sol(:,1);
sol=sol(:,2:5);


cheb_coeff = get_cheb_coeffs(sol, params);
x.a1=cheb_coeff(:,1)';
x.a2=cheb_coeff(:,2)';
x.a3=cheb_coeff(:,3)';
x.a4=cheb_coeff(:,4)';

% Perform Newton's method 

new_x = refine_cheb_orbit(x, mflds, params);
% new_x.a1 gives bs

phi_cheb = new_x.a1;
%%

phi = chebfun(1);
phi.domain = [-1,1];
phi.funs{1,1}.onefun.coeffs = [phi_cheb(1),2*phi_cheb(2:end)]';



for i = 1:2

% chebfun integrator

if i == 1
    v_u = real(evecs.u(:,1));
else
    v_u = imag(evecs.u(:,1));
end


N = chebop(-1,1);
N.op = @(t,h1,h2,h3,h4) [diff(h1)-L1*(h4);
                   diff(h2)-L1*(h3-2*h4);
                   diff(h3)-L1*(-h1+(2*params.nu*phi-3*phi^2-params.mu)*h1);
                   diff(h4)-L1*(h2)];

N.lbc = v_u;
[h1,h2,h3,h4] = N\0;


% newton

n = length(h1);
NN = n;
phi_cheb = [phi_cheb, zeros(1,NN)];
phi_cheb = phi_cheb(1:NN);

h = zeros(1,4*NN);
h(1:n) = chebcoeffs(h1)/2;
h(1) = h(1)*2;
h(NN+1:NN+n) = chebcoeffs(h2)/2;
h(NN+1) = h(NN+1)*2;
h(2*NN+1:2*NN+n) = chebcoeffs(h3)/2;
h(2*NN+1) = h(2*NN+1)*2;
h(3*NN+1:3*NN+n) = chebcoeffs(h4)/2;
h(3*NN+1) = h(3*NN+1)*2;


for j = 1:5

    h = h - (chebDF(phi_cheb,NN,params)\chebF(h,v_u,phi_cheb,NN,params))';

end

disp('norm of F(h) at end of Newton:')
disp(norm(chebF(h,v_u,phi_cheb,NN,params)))

if i == 1  %make one of these phi'                                                                                                                                                                      %hi
    H.a11 = chebcoeff_to_function(h(1:NN));
    H.a21 = chebcoeff_to_function(h(NN+1:2*NN));
    H.a31 = chebcoeff_to_function(h(2*NN+1:3*NN));
    H.a41 = chebcoeff_to_function(h(3*NN+1:4*NN));
else
    H.a12 = chebcoeff_to_function(h(1:NN));
    H.a22 = chebcoeff_to_function(h(NN+1:2*NN));
    H.a32 = chebcoeff_to_function(h(2*NN+1:3*NN));
    H.a42 = chebcoeff_to_function(h(3*NN+1:4*NN));
end


end



detA = (H.a11 .* H.a42) - (H.a12 .* H.a41);
plot(L1*(-1:.05:1),detA)





