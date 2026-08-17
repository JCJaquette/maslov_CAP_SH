function x = cauchyProd2D_mn(a,b,m,n)
% cauchy product for 2D Taylor series coeffs

    aS = size(a); bS = size(b);
    if isintval(a) || isintval(b)
        z = intval(1)*zeros(m+1,n+1);
    else
        z = zeros(m+1,n+1);
    end
    z(1:aS(1),1:aS(2)) = a;
    a = z;
    z = 0*z;
    z(1:bS(1),1:bS(2)) = b;
    b = z;%in case m or n is too big pad matrix a and b with zeros

    amat = a(m+1:-1:1,n+1:-1:1);
    bmat = b(1:m+1,1:n+1);

    x = sum(sum(amat.*bmat));

end

