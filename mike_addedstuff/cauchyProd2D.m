function ab = cauchyProd2D(a,b)


    aL = length(a(1,:)); aH = length(a(:,1));
    bL = length(b(1,:)); bH = length(b(:,1));   
    m = aH+bH;
    n = aL+bL;
    ab = zeros(m-1,n-1);

    for i = 1:m
        for j = 1:n

            ab(i,j) = cauchyProd2D_mn(a,b,i-1,j-1); 

        end
    end

end

