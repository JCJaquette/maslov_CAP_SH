function ab = cauchyProd2D(a,b)


    aL = length(a(1,:)); aH = length(a(:,1));
    bL = length(b(1,:)); bH = length(b(:,1));   
    m = aH+bH;
    n = aL+bL;
    
    if isintval(a) || isintval(b)
        ab = intval(1)*zeros(m-1,n-1);
    else
        ab = zeros(m-1,n-1);
    end

    for i = 1:m-1
        for j = 1:n-1

            ab(i,j) = cauchyProd2D_mn(a,b,i-1,j-1); 

        end
    end

end

