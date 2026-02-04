function [point] = get_manifold_point(coeff,phi1,phi2,order)

    p=60;
    
    %We just get the one point we want from the plot_manifold function
    r=linspace(0,1,p);
    theta=linspace(0,2*pi, p);

    ps1s2 = zeros(4,1);
    s1=phi1;
    s2=phi2;

    for n=0:order
        for m=0:n

            point=reshape(coeff(n-m+1,m+1,:),[4,1]);            
            ps1s2=ps1s2+point.*(s1+1i*s2)^(n-m).*(s1-1i*s2)^m;

        end       
    end 

    point=ps1s2;

    
end

