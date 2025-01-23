function [out] = chebcProd_nth(a,b,k)
%kth element of the convolution product a*b, k = 0 gives (a*b)_0, 1 gives (a*b)_1... 
    
    N = max(length(a),length(b));
    if N == length(a)
        a = [a;zeros(k+1,1)];
        b = [b;zeros(N+k+1 - length(b),1)];
    end
    if N == length(b)
        a = [a;zeros(N+k+1 - length(a),1)];
        b = [b;zeros(k+1,1)];
    end 

    out = dot(a(1:k+1),b(k+1:-1:1)) + dot(a(2:N+1),b(k+2:k+N+1)) + dot(b(2:N+1),a(k+2:k+N+1));

end