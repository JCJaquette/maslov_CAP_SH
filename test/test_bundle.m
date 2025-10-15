clear

tic 

x = zeros(1, 4); 
params.lambda = 0; 
params.mu = 0.05; 
params.nu = 1.6; 
params.scale = 1/7000;
params.isIntval = 0;

if params.isIntval
    zero=intval(0);
else
    zero=0;
end

order = 10-1; 
params.order = order; 
params.mfld.order = order; 



%  TODO: Make get_mflds universal
mflds=get_mflds(params);

BOOL_stable = 1;
[mflds,r_min_s]=mfld_poly(params, mflds,BOOL_stable );
 


% Computes Manifold coeff, and bundle Coeff.
[bndl] = getAllBundleCoefficients(params,mflds);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Compare 
% % % Define change of coordinates  symplectic = S * normal
    S= [1 0 0 0;
        0 0 1 0;
        0 2 0 1;
        0 1 0 1];
    S_inv...
     = [1 0 0 0;
        0 0 0 1;
        0 1 0 0;
        0 0 1 -2];

% Start with a point on the manifold 

psi_0=3.669582882882883;
    psi1 = cos(psi_0);
    psi2 = sin(psi_0); 

    sigma_0 = exp(1i*4);

    
    tmax= 30;
    tspan = [0,tmax];

   
    pointy = mfld_one_point(real(sigma_0),imag(sigma_0),mflds.stable.coeffs, params);
    pointy = real(pointy);
 
    
% Start with a tangent vectors

Vtilde0 = [1;0;0;0];

Bundle_0 = bndl_one_point(real(sigma_0), imag(sigma_0), bndl, params); %TODO

V0 = Bundle_0 *Vtilde0;

V0=V0;

V0_symplectic =S*V0;

V0_symplectic
W0=bigPhi(0,[0*V0;V0_symplectic],params)
W0(5:end)./V0_symplectic
% Something is going wrong here. The eigenvectors are not matching up with
% what they should be. Then we need to check that the unstable resonant
% bundles grow appropriately. 
return
% Integrate it with the normal form 

uV_0 = [sigma_0;conj(sigma_0);Vtilde0];

normal_form_ode(0,uV_0,params,mflds,bndl)
[t,uV] = ode45(@(t,uV)normal_form_ode(t,uV,params,mflds,bndl),tspan,uV_0);

V_list =uV;

for i = 1:length(t)
    phi_local = uV(i,1);
    bndl_local = bndl_one_point(real(phi_local), imag(phi_local), bndl, params);
    V_local = S*bndl_local*transpose(uV(i,3:6));
    V_list(i,3:6) = transpose(V_local );
end

% Integrate it with the variational equation

Y_0 = [pointy;V0_symplectic];


[t2,Y] = ode45(@(T,Y)bigPhi(T,Y,params),tspan,Y_0)


% % % % % Plot everything
    figure(2)
    % plot(t_list,point_list(:,2))
    
    bundle_no = 1;
    plot(t2,Y(:,4+bundle_no ))
    hold on 
    plot(t,V_list(:,2+bundle_no ))
    hold off
    legend('basic','Normal Form')
 



function duVdt = normal_form_ode(t,uV,params,mflds,bndl)
 

Omega = diag([mflds.values.s mflds.values.u]);
u=uV(1:2);

V = uV(3:6);

a_13_20 = bndl.normalForm(1,3,2+1,0+1);
a_14_11 = bndl.normalForm(1,3,1+1,1+1);
a_23_11 = bndl.normalForm(2,3,1+1,1+1);
a_24_02 = bndl.normalForm(2,4,0+1,0+1);



dudt = transpose(mflds.values.s) .* u;
dVdt = Omega*V + [(a_13_20*u(1)^2 + a_14_11*u(1)*u(2))*V(3);(a_23_11*u(1)*u(2)+a_24_02*u(2)^2)*V(4);0;0];
% TODO Add normal form
duVdt = [dudt;dVdt];

end
