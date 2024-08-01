function [sum] = dcp(a,b,c,k)

    lc = length(c);
    lb = length(b);
    sum = 0;
    a = [a;zeros(2*(lc+lb),1)];
    b = [b;zeros(2*(lc+lb),1)];
    c = [c;zeros(2*(lc+lb),1)];

    for n = -lc:lc
        if c(abs(n)+1) ~= 0
            for m = -lb:lb
    
                sum = sum + c(abs(n)+1)*b(abs(m)+1)*a(abs(k-n-m)+1);
    
            end
        end
    end

end

