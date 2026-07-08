function injective = check_A_injective(params,DF,Am)
    m = params.pulse.order;

    % Compute the (expensive) interval product I - A*DF only once.
    % Use the ell-1 operator norm (max absolute column sum): it is the
    % operator norm on the weighted ell-1 sequence space we work in, and it
    % is O(n^2) rather than the O(n^3) SVD that norm(X) (the 2-norm) needs.
    IminusADF_mag = mag(eye(4*m+3)-Am*DF);
    nrm = norm(IminusADF_mag,1);

    disp('Norm of I - ADF: ')
    disp(nrm)

    if nrm>=1
        disp('Stop! The matrix Am is not injective.');
        injective = 0;
        return
    else
        disp('Good to go! The matrix Am is injective.');
        injective=1;
    end
end
