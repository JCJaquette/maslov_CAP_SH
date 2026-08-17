function y = get_newton_seed(params,mflds,pulseplot)
    
    mflds = struct_intvaltodouble(mflds);
    params = struct_intvaltodouble(params);     

    psoln = getPulse(params);

    [u_pts,u_phi1phi2s] = get_mani_points(mflds.unstable.coeffs,params.mfld.order,1);
    s_pts = get_mani_points(mflds.stable.coeffs,params.mfld.order,1);

    manifold_dimensions = size(s_pts);

    thetas = linspace(0,2*pi,manifold_dimensions(1));
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
    righthalf = psoln(n+1:end,2:end); % Note, righthalf also includes the center

    % Start from the middle, find the first time it gets close enough to
    % the manifold
    for k_s = 1:n+1  
        [closest_pts_on_s(k_s,:), distances_s(k_s), manifold_index_s(k_s,:)] = closestpt(s_pts,(righthalf(k_s,:))');
    end
    for k_s = 1:n+1
        if distances_s(k_s) < boundary_distance 
            break
        end
    end

    for k_u = 1:n  
        [closest_pts_on_u(k_u,:), distances_u(k_u), manifold_index_u(k_u,:)] = closestpt(u_pts,(lefthalf(k_u,:))');
    end
    % Start from the middle, find the first time it gets close enough to
    % the manifold
    for k_u = n:-1:1
        if distances_u(k_u) < boundary_distance
            break
        end
    end

    L_mns = psoln(k_u,1);
    L_pls = psoln(n+k_s,1); % Also needs to account for the middle point

    % Define k as the offset from the middle point
    if abs(L_pls) >= abs(L_mns)
        L = abs(L_pls);
        k = k_s-1;
    else
        L = abs(L_mns);
        k = (n-k_u)+1;
    end

    k_half_ind_left = n-k+1;
    k_half_ind_right = k+1;

    % time_left = psoln(n+1-k)
    % time_right = psoln(n+1+k)

    Lsoln = psoln(n+1-k:n+1+k,2:end);

    phi1 = u_phi1phi2s(manifold_index_u(k_half_ind_left,1),manifold_index_u(k_half_ind_left,2),1);
    phi2 = u_phi1phi2s(manifold_index_u(k_half_ind_left,1),manifold_index_u(k_half_ind_left,2),2);

    y = chebfuncoeffs(Lsoln,params.pulse.order);
    y.phi1 = phi1;
    y.phi2 = phi2;
    y.psi = thetas(manifold_index_s(k_half_ind_right,2));
    y.Lbvp = L;


if pulseplot

%%%%%%%%%%%%%%%%%
% Generate Plots 
%%%%%%%%%%%%%%%%

y1 = chebcoeff_to_function(y.a1);
y2 = chebcoeff_to_function(y.a2);
y3 = chebcoeff_to_function(y.a3);
y4 = chebcoeff_to_function(y.a4);

dom = -1:.05:1;

raw_times = psoln(n+1-k:n+1+k,1);
raw_times_norm = (raw_times - raw_times(1)) / (raw_times(end) - raw_times(1)) * 2 - 1;

yo1 = interp1(raw_times_norm, Lsoln(:,1), dom, 'pchip');
yo2 = interp1(raw_times_norm, Lsoln(:,2), dom, 'pchip');
yo3 = interp1(raw_times_norm, Lsoln(:,3), dom, 'pchip');
yo4 = interp1(raw_times_norm, Lsoln(:,4), dom, 'pchip');

 figure
 tiledlayout(4,1)
 nexttile
 plot(dom,y1)
 nexttile
 plot(dom,y2)
 nexttile
 plot(dom,y3)
 nexttile
 plot(dom,y4)
 title('Solution obtained via Newtons method.')


figure 
hold on 
plot(dom, y1, linewidth = 1.5, color = "#A2142F")
plot(dom(1), y1(1), 'o', 'MarkerFaceColor', "#A2142F")
plot(dom(end), y1(end), 'o', 'MarkerFaceColor', "#A2142F")
hold off
%xlabel('$t$', Interpreter = 'latex', FontSize=14)
%ylabel('$\varphi(t)$', Interpreter = 'latex', FontSize=14)

 
 
 dom = L*dom;
 
 figure 
 tiledlayout(4,1)
 nexttile
 hold on
 plot(dom,yo1)
 plot(dom,y1)
 legend('Chebyshev Rep.', 'Refined Chebyshev Rep.')
 hold off
 nexttile
 hold on
 plot(dom,yo2)
 plot(dom,y2)
 hold off
 nexttile
 hold on
 plot(dom,yo3)
 plot(dom,y3)
 hold off
 nexttile
 hold on 
 plot(dom,yo4)
 plot(dom,y4)
 hold off 
 
 
figure
tiledlayout(2,1)
nexttile
plot_coeff(mflds.unstable.coeffs, params.mfld.order);
title('Unstable Mfld Coeff.')
nexttile
plot_coeff(mflds.stable.coeffs, params.mfld.order);
title('Stable Mfld Coeff.')


mflds.pts.s=mfld_points(mflds.stable.coeffs, params);
mflds.pts.u=mfld_points(mflds.unstable.coeffs, params);

figure
figure
hold on
surf(mflds.pts.s(:,:,1),mflds.pts.s(:,:,2),mflds.pts.s(:,:,4), 'FaceColor','r', 'FaceAlpha',0.5, 'EdgeColor','none');  
xlabel('x1');
ylabel('x2');
zlabel('x4');
surf(mflds.pts.u(:,:,1),mflds.pts.u(:,:,2),mflds.pts.u(:,:,4), 'FaceColor','g', 'FaceAlpha',0.5, 'EdgeColor','none');  
plot3(y1,y2,y4);
text(y1(1), y2(1), y4(1), '  P(\theta)', 'Color', 'g', 'FontSize', 12, 'Interpreter', 'tex');
text(y1(end), y2(end), y4(end), '  Q(\phi)', 'Color', 'r', 'FontSize', 12, 'Interpreter', 'tex');
title('Stable and Unstable Manifolds');
legend('Stable','Unstable','Hom. Orbit');
hold off



figure
grid on 
hold on
colormap(spring)
surf(mflds.pts.s(:,:,1),mflds.pts.s(:,:,2),mflds.pts.s(:,:,4), mflds.pts.s(:,:,3), 'FaceAlpha',0.4, 'EdgeColor', '#554b1c');
xlabel('$x_1$', Interpreter = 'latex');
ylabel('$x_2$', Interpreter = 'latex');
zlabel('$x_4$', Interpreter = 'latex')
surf(mflds.pts.u(:,:,1),mflds.pts.u(:,:,2),mflds.pts.u(:,:,4),'FaceAlpha',0.4);  
plot3(y1,y2,y4, lineWidth = 2, Color="#A2142F");
plot3(y1(1), y2(1), y4(1), 'o', 'MarkerFaceColor', "#A2142F")
plot3(y1(end), y2(end), y4(end), 'o', 'MarkerFaceColor', "#A2142F")
text(y1(1), y2(1), y4(1), '  P(\theta)', 'Color', '#7E2F8E', 'FontSize', 12, 'Interpreter', 'tex');
text(y1(end), y2(end), y4(end), '  Q(\phi)', 'Color', '#EDB120', 'FontSize', 12, 'Interpreter', 'tex');
%title('Validated Manifolds and $\varphi(x)$ Trajectory from ');
legend('Stable manifold','Unstable manifold','$\varphi(x)$', Interpreter = 'latex');
hold off

% color_stable = '#554b1c';
% color_unstable = "#7E2F8E";
color_stable = 'b';
color_unstable = 'r';
figure
hold on
surf(mflds.pts.s(:,:,1),mflds.pts.s(:,:,2),mflds.pts.s(:,:,4), 'FaceAlpha',0.4, 'FaceColor', color_stable , edgeColor = "none");
xlabel('$x_1$', Interpreter = 'latex');
ylabel('$x_2$', Interpreter = 'latex');
zlabel('$x_4$', Interpreter = 'latex')
surf(mflds.pts.u(:,:,1),mflds.pts.u(:,:,2),mflds.pts.u(:,:,4),'FaceAlpha',0.4, 'FaceColor',color_unstable, edgeColor = "none"); 
plot3(y1,y2,y4, lineWidth = 2, Color="#A2142F");
plot3(y1(1), y2(1), y4(1), 'o', 'MarkerFaceColor', "#A2142F")
plot3(y1(end), y2(end), y4(end), 'o', 'MarkerFaceColor', "#A2142F")
text(y1(end), y2(end), y4(end), '  P(\theta)', 'Color', color_stable, 'FontSize', 12, 'Interpreter', 'tex');
text(y1(1), y2(1), y4(1), '  Q(\phi)', 'Color', color_unstable, 'FontSize', 12, 'Interpreter', 'tex');
%title('Validated Manifolds and $\varphi(x)$ Trajectory from ');
legend('Stable manifold','Unstable manifold','$\varphi(x)$', Interpreter = 'latex');
hold off

end

end 