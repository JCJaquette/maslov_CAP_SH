clear
close all
%1 is mu=.05,nu=1.6, 0 branch, 2 is mu=.05,nu=1.6, π branch, 3 is mu=.2,nu=1.6

    n = 1;

    [params,mflds,psoln] = getparamsVerify(n);


    hold on
    [u_pts,u_phi1phi2s] = plot_manifold(mflds.unstable.coeffs,25,'red');
    s_pts = plot_manifold(mflds.stable.coeffs,25,'blue');
    thetas = linspace(0,2*pi,60);
    boundary_distance = avgnorms(u_pts);
    boundary_distance = boundary_distance*params.bd_scale;

   
    n = floor(length(psoln(:,2))/2);
    distances_u = zeros(n,1);
    closest_pts_on_u = zeros(n,4);
    distances_s = zeros(n,1);
    closest_pts_on_s = zeros(n,4);
    manifold_index_s = zeros(n,2);
    manifold_index_u = zeros(n,2);    

    
    lefthalf = psoln(1:n,2:end);
    righthalf = psoln(n+1:end,2:end);

    for k_s = 1:n  
        [closest_pts_on_s(k_s,:), distances_s(k_s), manifold_index_s(k_s,:)] = closestpt(s_pts,(righthalf(k_s,:))');
    end
    for k_s = 1:n
        if distances_s(k_s) < boundary_distance 
            break
        end
    end

    for k_u = 1:n  
        [closest_pts_on_u(k_u,:), distances_u(k_u), manifold_index_u(k_u,:)] = closestpt(u_pts,(lefthalf(k_u,:))');
    end
    for k_u = n:-1:1
        if distances_u(k_u) < boundary_distance
            break
        end
    end

    L_mns = psoln(k_u,1);
    L_pls = psoln(n+1+k_s,1);
    if abs(L_pls) > abs(L_mns)
        params.L = abs(L_pls);
        k = k_s;
    else
        params.L = abs(L_mns);
        k = n-k_u+1;
    end

    k_half_ind_left = n-k+1;
    k_half_ind_right = k;

    lefttime = psoln(n-k+1,1);
    righttime = psoln(n+k+1,1);

    Lsoln = psoln(n-k+1:n+k+1,2:end);

    left_endpt_u = closest_pts_on_u(k_half_ind_left,:);    
    right_endpt_s = closest_pts_on_s(k_half_ind_right,:);

    phi1 = u_phi1phi2s(manifold_index_u(k_half_ind_left,1),manifold_index_u(k_half_ind_left,2),1);
    phi2 = u_phi1phi2s(manifold_index_u(k_half_ind_left,1),manifold_index_u(k_half_ind_left,2),2);

    y = chebfuncoeffs(Lsoln,params.cheb.order);
    y.phi1 = phi1;
    y.phi2 = phi2;
    y.psi = thetas(manifold_index_s(k_half_ind_right,2));
    
    plot3(Lsoln(:,1),Lsoln(:,2),Lsoln(:,4),'black','LineWidth',1)
    
    plot3(right_endpt_s(1),right_endpt_s(2),right_endpt_s(4),'. black','MarkerSize',16);
    plot3(left_endpt_u(1),left_endpt_u(2),left_endpt_u(4),'. black','MarkerSize',16);
    
    plot3(Lsoln(1,1),Lsoln(1,2),Lsoln(1,4),'. black','MarkerSize',16)
    plot3(Lsoln(end,1),Lsoln(end,2),Lsoln(end,4),'. black','MarkerSize',16)
    

    new_y = refine_cheb_orbit(y,mflds,params);

    yo1 = chebcoeff_to_function(new_y.a1);
    yo2 = chebcoeff_to_function(new_y.a2);
    yo3 = chebcoeff_to_function(new_y.a3);
    yo4 = chebcoeff_to_function(new_y.a4);
    plot3(yo1 ,yo2 ,yo4 ,'LineWidth',1)

    legend('Unstable Manifold','Stable Manifold', 'Pulse')

    
    verify_homoclinic_orbit(params,mflds,new_y,params.new);


