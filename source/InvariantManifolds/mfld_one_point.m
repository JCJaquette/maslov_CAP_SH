function val = mfld_one_point(phi1, phi2, coeff, params)
    order=params.mfld.order;
    dim = length(coeff(1,1,:));
    val=zeros(4,1);
    for n=0:order
        for m=0:n
            point=reshape(coeff(n-m+1,m+1,:),[dim,1]);
            val=val+point.*(phi1+1i*phi2)^(n-m)*(phi1-1i*phi2)^m;
        end
    end
       if norm(val-real(val)) > 1e-8%10->8
            disp(imag(val))
            msg = 'Error occurred. The manifold is complex valued.';
            error(msg);
       end
end 

