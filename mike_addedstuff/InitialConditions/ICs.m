clear

%1 is mu=.05,nu=1.6, 0 branch, 2 is mu=.05,nu=1.6, π branch, 3 is mu=.2,nu=1.6
n = 2;



if n == 1

    params.rho = 1 - .01;
    params.scale = 2e-1;
    params.mu=0.05;
    params.nu=1.6;
    params.lambda = 0;
    params.cheb.order=450;
    params.mfld.order=25;
    params.tol=4e-16;
    
    load("psoln1.mat");

    mflds = get_mflds(params);
    hold on
    [u_pts,u_p1p2s] = plot_manifold(mflds.unstable.coeffs,25,'red');
    s_pts = plot_manifold(mflds.stable.coeffs,25,'blue');
    thetas = linspace(0,2*pi,60);
    boundary_distance = avgnorms(u_pts);
    boundary_distance = boundary_distance*.9;
    
    n = length(pSoln1.normalForm.sol(:,1));
    d_u = zeros(n,1);
    pts_u = zeros(n,4);
    d_s = zeros(n,1);
    pts_s = zeros(n,4);
    mps = zeros(n,2);
    
    for k = 1:n  
        [pts_s(k,:), d_s(k), mps(k,:)] = closestpt(s_pts,(pSoln1.normalForm.sol(k,:))');
    end
    for k = 1:n
        if d_s(k) < boundary_distance 
            break
        end
    end
    
    L_pl = pSoln1.normalForm.time(k);
    params.L = L_pl;
    pt_s = pts_s(k,:);
    Lsoln = pSoln1.normalForm.sol(1:k,:);    

    pSoln1.normalForm.sol(:,2) = -pSoln1.normalForm.sol(:,2);
    pSoln1.normalForm.sol(:,4) = -pSoln1.normalForm.sol(:,4);

    for k = 1:n  
        [pts_u(k,:), d_u(k)] = closestpt(u_pts,(pSoln1.normalForm.sol(k,:))');
    end
    for k = 1:n
        if d_u(k) < boundary_distance
            break
        end
    end

    L_mn = -pSoln1.normalForm.time(k);
    lefthalf = Lsoln(end:-1:1,:);

    righthalf = pSoln1.normalForm.sol(2:k,:);
    righthalf(:,2) = -righthalf(:,2);
    righthalf(:,4) = -righthalf(:,4);
    % Lsoln = [ lefthalf;
    %           righthalf];

    pt_u = pts_u(k,:);

    % phis
    phi1 = u_p1p2s(mps(k,1),mps(k,2),1);
    phi2 = u_p1p2s(mps(k,1),mps(k,2),2);

    y = chebfuncoeffs(Lsoln);
    y.phi1 = phi1;
    y.phi2 = phi2;
    y.psi = thetas(mps(k,2));
    
    plot3(Lsoln(:,1),Lsoln(:,2),Lsoln(:,4),'black','LineWidth',1)
    
    plot3(pt_s(1),pt_s(2),pt_s(4),'. black','MarkerSize',16);
    plot3(pt_u(1),pt_u(2),pt_u(4),'. black','MarkerSize',16);
    
    plot3(Lsoln(1,1),Lsoln(1,2),Lsoln(1,4),'. black','MarkerSize',16)
    plot3(Lsoln(end,1),Lsoln(end,2),Lsoln(end,4),'. black','MarkerSize',16)

    legend('Unstable Manifold','Stable Manifold', 'Pulse')
    

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if n == 2

    params.rho = 1 - .01;
    params.scale = 3e-1;
    params.mu=0.05;
    params.nu=1.6;
    params.lambda = 0;
    params.cheb.order=450;
    params.mfld.order=25;
    params.tol=4e-16;
    
    load("psoln2.mat");

    mflds = get_mflds(params);
    hold on
    [u_pts,u_p1p2s] = plot_manifold(mflds.unstable.coeffs,25,'red');
    s_pts = plot_manifold(mflds.stable.coeffs,25,'blue');
    thetas = linspace(0,2*pi,60);

    boundary_distance = avgnorms(u_pts);
    boundary_distance = boundary_distance*.7;
    
    n = length(pSoln2.normalForm.sol(:,1));
    d_u = zeros(n,1);
    pts_u = zeros(n,4);
    d_s = zeros(n,1);
    pts_s = zeros(n,4);
    mps = zeros(n,2);
    
    for k = 1:n  
        [pts_s(k,:), d_s(k), mps(k,:)] = closestpt(s_pts,(pSoln2.normalForm.sol(k,:))');
    end
    for k = 1:n
        if d_s(k) < boundary_distance 
            break
        end
    end
 
    point_on_manifold = pts_s(k,:);
    point_on_pulse = pSoln2.normalForm.sol(k,:);

    L_pl = pSoln2.normalForm.time(k);
    params.L = L_pl;
    pt_s = pts_s(k,:);
    Lsoln = pSoln2.normalForm.sol(1:k,:);    
    
    pSoln2.normalForm.sol(:,2) = -pSoln2.normalForm.sol(:,2);
    pSoln2.normalForm.sol(:,4) = -pSoln2.normalForm.sol(:,4);

    for k = 1:n  
        [pts_u(k,:), d_u(k)] = closestpt(u_pts,(pSoln2.normalForm.sol(k,:))');
    end
    for k = 1:n
        if d_u(k) < boundary_distance
            break
        end
    end

    L_mn = -pSoln2.normalForm.time(k);
    lefthalf = Lsoln(end:-1:1,:);

    righthalf = pSoln2.normalForm.sol(2:k,:);
    righthalf(:,2) = -righthalf(:,2);
    righthalf(:,4) = -righthalf(:,4);
    Lsoln = [ lefthalf;
             righthalf];
    pt_u = pts_u(k,:);

    % phis
    phi1 = u_p1p2s(mps(k,1),mps(k,2),1);
    phi2 = u_p1p2s(mps(k,1),mps(k,2),2);

    pSoln2.normalForm.sol(:,2) = -pSoln2.normalForm.sol(:,2);
    pSoln2.normalForm.sol(:,4) = -pSoln2.normalForm.sol(:,4);

    y = chebfuncoeffs(Lsoln);
    y.phi1 = phi1;
    y.phi2 = phi2;
    y.psi = thetas(mps(k,2));

    plot3(Lsoln(:,1),Lsoln(:,2),Lsoln(:,4))
    
    plot3(pt_s(1),pt_s(2),pt_s(4),'. black','MarkerSize',16);
    plot3(pt_u(1),pt_u(2),pt_u(4),'. black','MarkerSize',16);
    
    plot3(Lsoln(1,1),Lsoln(1,2),Lsoln(1,4),'. black','MarkerSize',16)
    plot3(Lsoln(end,1),Lsoln(end,2),Lsoln(end,4),'. black','MarkerSize',16)
    
    % figure 
    % hold on
    % plot(chebcoeff_to_function(y.a1))
    % plot(Lsoln(:,1))
    % 
    % figure 
    % hold on
    % plot(chebcoeff_to_function(y.a2))
    % plot(Lsoln(:,2))


end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if n == 3

    params.rho = 1 - .01;
    params.scale = .9;
    params.mu=0.2;
    params.nu=1.6;
    params.lambda = 0;
    params.cheb.order=450;
    params.mfld.order=25;
    params.tol=4e-15;
    
    load("psoln3.mat");

    mflds = get_mflds(params);
    hold on
    [u_pts,u_p1p2s] = plot_manifold(mflds.unstable.coeffs,25,'red');
    s_pts = plot_manifold(mflds.stable.coeffs,25,'blue');
    thetas = linspace(0,2*pi,60);
    
    boundary_distance = avgnorms(u_pts);
    boundary_distance = boundary_distance*.9;
    
    n = length(pSoln3.normalForm.sol(:,1));
    d_u = zeros(n,1);
    pts_u = zeros(n,4);
    d_s = zeros(n,1);
    pts_s = zeros(n,4);
    mps = zeros(n,2);
    
    for k = 1:n  
        [pts_s(k,:), d_s(k), mps(k,:)] = closestpt(s_pts,(pSoln3.normalForm.sol(k,:))');
    end
    for k = 1:n
        if d_s(k) < boundary_distance 
            break
        end
    end
    
    L_pl = pSoln3.normalForm.time(k);
    params.L = L_pl;
    pt_s = pts_s(k,:);
    Lsoln = pSoln3.normalForm.sol(1:k,:);    
    
    pSoln3.normalForm.sol(:,2) = -pSoln3.normalForm.sol(:,2);
    pSoln3.normalForm.sol(:,4) = -pSoln3.normalForm.sol(:,4);

    for k = 1:n  
        [pts_u(k,:), d_u(k)] = closestpt(u_pts,(pSoln3.normalForm.sol(k,:))');
    end
    for k = 1:n
        if d_u(k) < boundary_distance
            break
        end
    end

    L_mn = -pSoln3.normalForm.time(k);
    Lsoln = [           Lsoln(end:-1:1,:);
             pSoln3.normalForm.sol(2:k,:)];
    pt_u = pts_u(k,:);

    % phis
    phi1 = u_p1p2s(mps(k,1),mps(k,2),1);
    phi2 = u_p1p2s(mps(k,1),mps(k,2),2);

    pSoln3.normalForm.sol(:,2) = -pSoln3.normalForm.sol(:,2);
    pSoln3.normalForm.sol(:,4) = -pSoln3.normalForm.sol(:,4);

    y = chebfuncoeffs(Lsoln);
    y.phi1 = phi1;
    y.phi2 = phi2;
    y.psi = thetas(mps(k,2));
    
    plot3(Lsoln(:,1),Lsoln(:,2),Lsoln(:,4))
    
    plot3(pt_s(1),pt_s(2),pt_s(4),'. black','MarkerSize',16);
    plot3(pt_u(1),pt_u(2),pt_u(4),'. black','MarkerSize',16);
    
    plot3(Lsoln(1,1),Lsoln(1,2),Lsoln(1,4),'. black','MarkerSize',16)
    plot3(Lsoln(end,1),Lsoln(end,2),Lsoln(end,4),'. black','MarkerSize',16)

end
