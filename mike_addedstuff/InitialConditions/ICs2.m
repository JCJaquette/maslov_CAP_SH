clear

%1 is mu=.05,nu=1.6, 0 branch, 2 is mu=.05,nu=1.6, π branch, 3 is mu=.2,nu=1.6
n = 1;



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
    manifold_index_s = zeros(n,2);
    manifold_index_u = zeros(n,2);    

    
    for k1 = 1:n  
        [pts_s(k1,:), d_s(k1), manifold_index_s(k1,:)] = closestpt(s_pts,(pSoln1.normalForm.sol(k1,:))');
    end
    for k1 = 1:n
        if d_s(k1) < boundary_distance 
            break
        end
    end

    pSoln1.normalForm.sol(:,2) = -pSoln1.normalForm.sol(:,2);
    pSoln1.normalForm.sol(:,4) = -pSoln1.normalForm.sol(:,4);

    for k2 = 1:n  
        [pts_u(k2,:), d_u(k2), manifold_index_u(k2,:)] = closestpt(u_pts,(pSoln1.normalForm.sol(k2,:))');
    end
    for k2 = 1:n
        if d_u(k2) < boundary_distance
            break
        end
    end
    pSoln1.normalForm.sol(:,2) = -pSoln1.normalForm.sol(:,2);
    pSoln1.normalForm.sol(:,4) = -pSoln1.normalForm.sol(:,4);

    k = min(k1,k2);

    params.L = pSoln1.normalForm.time(k);
    pt_u = pts_u(k,:);    
    pt_s = pts_s(k,:);

    
    lefthalf = pSoln1.normalForm.sol(1:k,:);   
    lefthalf = lefthalf(end:-1:1,:);
    lefthalf(:,2) = -lefthalf(:,2);
    lefthalf(:,4) = -lefthalf(:,4);

    righthalf = pSoln1.normalForm.sol(2:k,:);

    Lsoln = [ lefthalf;
              righthalf];

    phi1 = u_p1p2s(manifold_index_u(k,1),manifold_index_u(k,2),1);
    phi2 = u_p1p2s(manifold_index_u(k,1),manifold_index_u(k,2),2);

    y = chebfuncoeffs(Lsoln);
    y.phi1 = phi1;
    y.phi2 = phi2;
    y.psi = thetas(manifold_index_s(k,2));
    
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
    thetas = linspace(0,2*pi,100);
    boundary_distance = avgnorms(u_pts);
    boundary_distance = boundary_distance*.9;
    
    n = length(pSoln2.normalForm.sol(:,1));
    d_u = zeros(n,1);
    pts_u = zeros(n,4);
    d_s = zeros(n,1);
    pts_s = zeros(n,4);
    manifold_index_s = zeros(n,2);
    
    for k1 = 1:n  
        [pts_s(k1,:), d_s(k1), manifold_index_s(k1,:)] = closestpt(s_pts,(pSoln2.normalForm.sol(k1,:))');
    end
    for k1 = 1:n
        if d_s(k1) < boundary_distance 
            break
        end
    end

    pSoln2.normalForm.sol(:,2) = -pSoln2.normalForm.sol(:,2);
    pSoln2.normalForm.sol(:,4) = -pSoln2.normalForm.sol(:,4);

    for k2 = 1:n  
        [pts_u(k2,:), d_u(k2)] = closestpt(u_pts,(pSoln2.normalForm.sol(k2,:))');
    end
    for k2 = 1:n
        if d_u(k2) < boundary_distance
            break
        end
    end
    pSoln2.normalForm.sol(:,2) = -pSoln2.normalForm.sol(:,2);
    pSoln2.normalForm.sol(:,4) = -pSoln2.normalForm.sol(:,4);

    k = max(k1,k2);

    params.L = pSoln2.normalForm.time(k);
    pt_u = pts_u(k,:);    
    pt_s = pts_s(k,:);

    
    lefthalf = pSoln2.normalForm.sol(1:k,:);   
    lefthalf = lefthalf(end:-1:1,:);
    lefthalf(:,2) = -lefthalf(:,2);
    lefthalf(:,4) = -lefthalf(:,4);

    righthalf = pSoln2.normalForm.sol(2:k,:);

    Lsoln = [ lefthalf;
              righthalf];

    phi1 = u_p1p2s(manifold_index_s(k,1),manifold_index_s(k,2),1);
    phi2 = u_p1p2s(manifold_index_s(k,1),manifold_index_s(k,2),2);

    y = chebfuncoeffs(Lsoln);
    y.phi1 = phi1;
    y.phi2 = phi2;
    y.psi = thetas(manifold_index_s(k,2));
    
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
    manifold_index_s = zeros(n,2);
    
    for k1 = 1:n  
        [pts_s(k1,:), d_s(k1), manifold_index_s(k1,:)] = closestpt(s_pts,(pSoln3.normalForm.sol(k1,:))');
    end
    for k1 = 1:n
        if d_s(k1) < boundary_distance 
            break
        end
    end

    pSoln3.normalForm.sol(:,2) = -pSoln3.normalForm.sol(:,2);
    pSoln3.normalForm.sol(:,4) = -pSoln3.normalForm.sol(:,4);

    for k2 = 1:n  
        [pts_u(k2,:), d_u(k2)] = closestpt(u_pts,(pSoln3.normalForm.sol(k2,:))');
    end
    for k2 = 1:n
        if d_u(k2) < boundary_distance
            break
        end
    end
    pSoln3.normalForm.sol(:,2) = -pSoln3.normalForm.sol(:,2);
    pSoln3.normalForm.sol(:,4) = -pSoln3.normalForm.sol(:,4);

    k = min(k1,k2);

    params.L = pSoln3.normalForm.time(k);
    pt_u = pts_u(k,:);    
    pt_s = pts_s(k,:);

    
    lefthalf = pSoln3.normalForm.sol(1:k,:);   
    lefthalf = lefthalf(end:-1:1,:);
    lefthalf(:,2) = -lefthalf(:,2);
    lefthalf(:,4) = -lefthalf(:,4);

    righthalf = pSoln3.normalForm.sol(2:k,:);

    Lsoln = [ lefthalf;
              righthalf];

    phi1 = u_p1p2s(manifold_index_s(k,1),manifold_index_s(k,2),1);
    phi2 = u_p1p2s(manifold_index_s(k,1),manifold_index_s(k,2),2);

    y = chebfuncoeffs(Lsoln);
    y.phi1 = phi1;
    y.phi2 = phi2;
    y.psi = thetas(manifold_index_s(k,2));
    
    plot3(Lsoln(:,1),Lsoln(:,2),Lsoln(:,4),'black','LineWidth',1)
    
    plot3(pt_s(1),pt_s(2),pt_s(4),'. black','MarkerSize',16);
    plot3(pt_u(1),pt_u(2),pt_u(4),'. black','MarkerSize',16);
    
    plot3(Lsoln(1,1),Lsoln(1,2),Lsoln(1,4),'. black','MarkerSize',16)
    plot3(Lsoln(end,1),Lsoln(end,2),Lsoln(end,4),'. black','MarkerSize',16)

    legend('Unstable Manifold','Stable Manifold', 'Pulse')
end
