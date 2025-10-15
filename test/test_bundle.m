clear

tic 

x = zeros(1, 4); 
params.lambda = 0; 
params.mu = 0.05; 
params.nu = 1.6; 
params.scale = 1/8;
params.isIntval = 0;

if params.isIntval
    zero=intval(0);
else
    zero=0;
end

order = 40-1; 
params.order = order; 
params.mfld.order = order; 

% Summary:
% Everything is good!!


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
        0 1 0 0];
    S_inv...
     = [1 0 0 0;
        0 0 0 1;
        0 1 0 0;
        0 0 1 -2];


 
% Start with a point on the manifold 

psi_0=3.669582882882883;
    psi1 = cos(psi_0);
    psi2 = sin(psi_0); 

    sigma_0 = exp(1i*2)*.1;

    
    tmax= 50;
    tspan = [0,tmax];

   
    pointy = mfld_one_point(real(sigma_0),imag(sigma_0),mflds.stable.coeffs, params);
    pointy = real(pointy);
 
    
% Start with a tangent vectors
a00=0;
a33=1;
Vtilde0 = [a00;conj(a00);a33;conj(a33)];

Bundle_0 = bndl_one_point(real(sigma_0), imag(sigma_0), bndl, params); %TODO

V0 = Bundle_0 *Vtilde0;

V0=V0;

V0_symplectic =S*V0;

% V0_symplectic
% W0=bigPhi(0,[0*V0;V0_symplectic],params)
% W0(5:end)./V0_symplectic
% Something is going wrong here. The eigenvectors are not matching up with
% what they should be. Then we need to check that the unstable resonant
% bundles grow appropriately. 
% return
% Integrate it with the normal form 

uV_0 = [sigma_0;conj(sigma_0);Vtilde0];

opts = odeset('RelTol',1e-8,'AbsTol',1e-10);

normal_form_ode(0,uV_0,params,mflds,bndl)
[t,uV] = ode45(@(t,uV)normal_form_ode(t,uV,params,mflds,bndl),tspan,uV_0 ,opts );

V_list =0*[uV(:,1:2) uV];

for i = 1:length(t)
    phi_local = uV(i,1);
    pt_local = mfld_one_point(real(phi_local),imag(phi_local),mflds.stable.coeffs, params);
    V_list(i,1:4)=pt_local ;

    bndl_local = bndl_one_point(real(phi_local), imag(phi_local), bndl, params);
    % V_local = S*bndl_local*transpose(uV(i,3:6));
    V_local = transpose(uV(i,3:6));
    V_list(i,5:8) = transpose(V_local );
end

% Integrate it with the variational equation

Y_0 = [pointy;V0_symplectic];


[t1,Y] = ode45(@(T,Y)bigPhi(T,Y,params),t,Y_0 ,opts );

Y_bundle = Y;
for i = 1:length(t1)
    phi_local = uV(i,1);
    bndl_local = bndl_one_point(real(phi_local), imag(phi_local), bndl, params);
    V_local = bndl_local\S_inv*transpose(Y_bundle(i,5:8));
    Y_bundle(i,5:8) = transpose(V_local );
end

% E=[mflds.vectors.s mflds.vectors.u];
% traj_1 = transpose(E\S_inv*transpose(Y(:,5:8 )));
% traj_2 = transpose(E\S_inv*transpose(V_list(:,3:6 )));


E=[mflds.vectors.s mflds.vectors.u];
traj_1 = Y_bundle(:,5:8);
traj_2 = V_list(:,5:8 );

% % % % % Plot everything
figure(1)
plot(t1,real(Y_bundle(:,1)))
    hold on 
    plot(t,real(V_list(:,1)))
hold off 

    figure(2)
    % plot(t_list,point_list(:,2))
    
tiledlayout(2,2) 
for  bundle_no = 1:4
    nexttile
    plot(t1,real(traj_1(:,bundle_no)))
    hold on 
    plot(t,real(traj_2(:,bundle_no)))
    hold off


    legend('basic','Normal Form')
end
 



function duVdt = normal_form_ode(t,uV,params,mflds,bndl)
 

Omega = diag([mflds.values.s mflds.values.u]);
u=uV(1:2);

V = uV(3:6);

a_13_20 = bndl.normalForm(1,3,2+1,0+1);
a_14_11 = bndl.normalForm(1,4,1+1,1+1);
a_23_11 = bndl.normalForm(2,3,1+1,1+1);
a_24_02 = bndl.normalForm(2,4,0+1,2+1);



dudt = transpose(mflds.values.s) .* u;
dVdt = Omega*V +  [a_13_20*u(1)^2*V(3) + a_14_11*u(1)*u(2)*V(4);a_23_11*u(1)*u(2)*V(3)+a_24_02*u(2)^2*V(4);0;0];
% TODO Add normal form
duVdt = [dudt;dVdt];

end
