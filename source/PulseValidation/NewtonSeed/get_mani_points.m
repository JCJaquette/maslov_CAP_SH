function [points,s1s2s]=get_mani_points(coeff,order,scale)
    p=60;

    r=linspace(0,scale,p);
    theta=linspace(0,2*pi, p);

    points=zeros(p,p,4);
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
        points(j,k,:)=real(ps1s2);
        end
    end

end