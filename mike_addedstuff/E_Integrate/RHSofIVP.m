function [phiprime] = RHSofIVP(pulse,phi,mu,nu)

    phiprime = 0*pulse;

    phiprime(1,:) = pulse(4,:);
    phiprime(2,:) = pulse(3,:) - 2*pulse(4,:);
    phiprime(3,:) = -pulse(1,:) + (2*nu*phi(1,:) - 3*phi(1,:).*phi(1,:) -mu) .* pulse(1,:);
    phiprime(4,:) = pulse(2,:);

    phiprime = L*phiprime;
end

