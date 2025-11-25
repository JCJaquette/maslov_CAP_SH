function [phiprime] = RHSofODE(pulse,mu,nu)
%eqn 9.4 in thesis

    phiprime = 0*pulse;

    % now RHS of ODE

    phiprime(1,:) = pulse(2,:);
    phiprime(2,:) = pulse(3,:);
    phiprime(3,:) = pulse(4,:);
    phiprime(4,:) = -2*pulse(3,:) - (mu + 1)*pulse(1,:) + nu*pulse(1,:).^2 - pulse(1,:).^3;

end

