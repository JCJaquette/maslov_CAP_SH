function A = IntegrateCheb(a)

    N = length(a);

    A = zeros(N+1,1);
    n = (1./(1:N))';

    a = [a;0;0];

    A(2:N+1) = n.*(a(1:N) - a(3:N+2))/2;


end

