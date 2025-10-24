clear all 

params.mu = 0.05; 
params.nu = 1.6; 
params.scale = 1/7;
params.lambda = 0;

order = 10-1;  % Manifold

params.order = order; % TODO Why do we have two orders?
params.mfld.order = order; 

% Add This boolean for intervals
params.isIntval = 0;



% -----------------------------------------------------------------------
% Now we calculate the coefficients of the parameterization for the stable
% and unstable manifold up to a desired order. 
 
% TODO: Replace convolution with Taylor_Convolution. ??

%  TODO: Make get_mflds universal
mflds=get_mflds(params);

% Pick point on manifold
psi_0=3.669582882882883;
    psi1 = cos(psi_0);
    psi2 = sin(psi_0); 

    sigma_1 = exp(1i*2);

    dt = .01;
    tmax= 20;
    t_list = 0:dt:tmax;
    tspan = [0,tmax];

    N_t = length(t_list);
    sigma_list = sigma_1 * exp(mflds.values.s(1) * t_list);

    point_list = zeros(N_t,4);
    for i =1:N_t
        sigma=sigma_list(i);
        pointy = mfld_one_point(real(sigma),imag(sigma),mflds.stable.coeffs, params);
        point_list(i,:)=real(pointy );
    end

[t,y] = ode45(@(t,x)SH_comp(t,x,params),tspan,point_list(1,:))

    figure(2)
    plot(t_list,point_list(:,2))
    hold on 
    plot(t,y(:,2))
hold off


    %%%%%%%%%%%%%%%%%%%%%%%%%  Unstable
     
    dt = .01;
    tmax= 20;
    t_list = 0:dt:tmax;
    tspan = [0,tmax];

    sigma_1 = exp(tmax*mflds.values.s(1));


   N_t = length(t_list);
    sigma_list = sigma_1 * exp(mflds.values.u(1) * t_list);

    point_list = zeros(N_t,4);
    for i =1:N_t
        sigma=sigma_list(i);
        pointy = mfld_one_point(real(sigma),imag(sigma),mflds.unstable.coeffs, params);
        point_list(i,:)=real(pointy );
    end

[t,y] = ode45(@(t,x)SH_comp(t,x,params),tspan,point_list(1,:))

    figure(3)
    plot(t_list,point_list(:,2))
    hold on 
    plot(t,y(:,2))
hold off
