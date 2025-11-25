function [Fout] = chebF(h,hu,b,N,params)
% makes the psi given in thesis pg231, eqn12.6
% zero of this gives h coeffs
% b is soln cheb coeffs, hu is h(-1), N is order of series

    h1 = [h(1:N)];
    h2 = [h(N+1:2*N)];
    h3 = [h(2*N+1:3*N)];
    h4 = [h(3*N+1:4*N)];


    Fai0 = zeros(4,1);
    alt = ones(1,N-1);
    for k = 1:N-1

        if mod(k,2) == 1
            alt(k) = -1;
        end

    end    

    Fai0(1) = h1(1) + 2*dot(alt,h1(2:N)) - hu(1);
    Fai0(2) = h2(1) + 2*dot(alt,h2(2:N)) - hu(2);
    Fai0(3) = h3(1) + 2*dot(alt,h3(2:N)) - hu(3);
    Fai0(4) = h4(1) + 2*dot(alt,h4(2:N)) - hu(4);

    c = make_c(h1,h2,h3,h4,b,params,N);    

    Fa1 = zeros(N-1,1);
    Fa2 = Fa1;
    Fa3 = Fa1;
    Fa4 = Fa1;

    for k = 1:N-1

        if k+2 <= N
            Fa1(k) = 2*k*h1(k+1) - params.L * (c(1,k) - c(1,k+2));
            Fa2(k) = 2*k*h2(k+1) - params.L * (c(2,k) - c(2,k+2));
            Fa3(k) = 2*k*h3(k+1) - params.L * (c(3,k) - c(3,k+2));
            Fa4(k) = 2*k*h4(k+1) - params.L * (c(4,k) - c(4,k+2));
        end
        if k+1 == N
            Fa1(k) = 2*k*h1(k+1) - params.L * c(1,k);
            Fa2(k) = 2*k*h2(k+1) - params.L * c(2,k);
            Fa3(k) = 2*k*h3(k+1) - params.L * c(3,k);
            Fa4(k) = 2*k*h4(k+1) - params.L * c(4,k);    
        end

    end

    Fout = [Fai0(1);Fa1;Fai0(2);Fa2;Fai0(3);Fa3;Fai0(4);Fa4];

end

