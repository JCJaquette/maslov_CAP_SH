function [sum] = chebcPProd_kth(a,b,c,k)
%kth element of the Cauchy product a*b*c, k = 0 gives (a*b*c)_0, 1 gives (a*b*c)_1...

    lc = length(c);
    lb = length(b);
    sum = 0;
    a = [a;zeros(2*(lc+lb),1)];
    b = [b;zeros(2*(lc+lb),1)];
    c = [c;zeros(2*(lc+lb),1)];

    for n = -lc:lc
        for m = -lb:lb

            sum = sum + c(abs(n)+1)*b(abs(m)+1)*a(abs(k-n-m)+1);

        end
    end

end

