function [plotpoints,s1s2s]=plot_manifold(coeff,order,color)
    p=60;
    
    % We plot them in polar coordinates 
    r=linspace(0,1,p);
    theta=linspace(0,2*pi, p);

    plotpoints=zeros(p,p,4);
    s1s2s = zeros(p,p,2);
    for j=1:p
        for k=1:p
            ps1s2 = zeros(4,1);
            s1=r(j)*cos(theta(k));
            s2=r(j)*sin(theta(k));
            s1s2s(j,k,:) = [s1;s2];
        for n=0:order
            for m=0:n
                point=reshape(coeff(n-m+1,m+1,:),[4,1]);
                
                ps1s2=ps1s2+point.*(s1+1i*s2)^(n-m).*(s1-1i*s2)^m;
            end
           
        end 
        plotpoints(j,k,:)=real(ps1s2);
        end
    end
    surf(plotpoints(:,:,1),plotpoints(:,:,2),plotpoints(:,:,4), 'FaceColor',color, 'FaceAlpha',0.5, 'EdgeColor','none');  

end