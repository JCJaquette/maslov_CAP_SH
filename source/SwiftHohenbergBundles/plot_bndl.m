function plots=plot_bndl(params,mflds,bndl,color)
% This function evaluates the function P on a grid of points and plots them
% Inputs: coeff - corresponds to the coefficients computed for either the
%                 stable or unstable manifold
%         order - order of the parameterization 
% Outputs: plots = 0
%          also generates a 3D plot. In this case, we omit the third
%          component (seemed to generate the best plot)
    

    All_Bundle_coeffs = bndl.coeffs;
    coeff = mflds.stable.coeffs;
    if params.isIntval
        coeff=coeff.mid;
        All_Bundle_coeffs=All_Bundle_coeffs.mid;
    end
    order = params.mfld.order; 

    % TODO: Add intlab compatibility
    
    % We plot them in polar coordinates 
    
    p=75;
    p_r= ceil(p/(2*pi));
    p_theta = p;

    r=linspace(0,1,p_r ).^.9;
    theta=linspace(0,2*pi, p_theta);

    plotpoints=zeros(p_r,p_theta,4);

    
    E = [ real(mflds.vectors.s(:,1)) imag(mflds.vectors.s(:,1)) real(mflds.vectors.u(:,1)) imag(mflds.vectors.u(:,1))];
    if params.isIntval
        E=E.mid;
    end
    % E=eye(4);
    bndle_no=3;
    for j=1:p_r
        for k=1:p_theta
            % pt_local = zeros(4,1);
            s1=r(j)*cos(theta(k));
            s2=r(j)*sin(theta(k));
            pt_local = E\mfld_one_point(s1,s2,coeff, params);
            % bndl_local = bndl_one_point(s1, s2, bndl, params);
             
            plotpoints(j,k,:)=real(pt_local);
            % plotbndl(j,k,:)=E\real(bndl_local(:,bndle_no));
            % plotbndl_norm(j,k)=norm(inv(bndl_local),1);
        end
    end
    X = plotpoints(:,:,1);
    Y = plotpoints(:,:,2);
    Z = plotpoints(:,:,4);

    p=40;
    p_r= ceil(p/(2*pi));
    p_theta = p;

    p_r2 = ceil(p_r/2);
    p_theta2 = ceil(p_theta/4);
    
    plotpoints2=zeros(p_r2,p_theta2,4);

    r2=linspace(1/(p_r2*2),1,p_r2 ).^.9;
    theta2=linspace(0,2*pi, p_theta2);

    plotbndl=zeros(p_r2,p_theta2,4);
    plotbndl_norm=zeros(p_r2,p_theta2);
    for j=1:p_r2
        for k=1:p_theta2
            pt_local = zeros(4,1);
            s1=r2(j)*cos(theta2(k));
            s2=r2(j)*sin(theta2(k));
            pt_local = E\mfld_one_point(s1,s2,coeff, params);
            bndl_local = bndl_one_point(s1, s2, bndl, params);

            if params.isIntval
                bndl_local =bndl_local .mid;
            end
             
            plotpoints2(j,k,:)=real(pt_local);
            plotbndl(j,k,:)=E\real(bndl_local(:,bndle_no));
            plotbndl_norm(j,k)=sum(abs((plotbndl(j,k,:))));
        end
    end

    X2 = plotpoints2(:,:,1);
    Y2 = plotpoints2(:,:,2);
    Z2 = plotpoints2(:,:,4);
    U = -plotbndl(:,:,1);
    V = -plotbndl(:,:,2);
    W = -plotbndl(:,:,4);


    q=quiver3(X2,Y2,Z2,U,V,W,'filled','r','LineWidth',.9,'MaxHeadSize',.06);
    
    hold on 
    % surf(X,Y,Z, 'FaceColor',color, 'FaceAlpha',.75,'EdgeColor','none');  
    surf(X,Y,Z, 'FaceColor',color, 'EdgeColor','none');  
    zlim([-.5,.125])
    camlight
    hold off 
    xlabel('$e^s_{Re}$','Interpreter','latex')
    ylabel('$e^s_{Im}$','Interpreter','latex')
    zlabel('$e^u_{Re}$','Interpreter','latex','Rotation',0)
    plots=0;
end