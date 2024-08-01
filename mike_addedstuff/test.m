clear
L=25;
tmin = -L;
tmax =L;
N = chebop(tmin,tmax);
nu = 1.6;
mu = .1;


A = [0,0,0,1;
    0,0,1,-2;
    -1-mu,0,0,0;
    0,1,0,0];
[V,D]=eig(A);
[~,ind]=sort(diag(D));
% Unstable eigenvector:
v_u = V(:,ind(end));
v_u=v_u/norm(v_u);

max_eig_real = max(real(eig(A)))
max_eig_real=0*max_eig_real;
phi =   chebfun(@(t) sech(t)*cos(t),[tmin,tmax]);
N.op = @(t,h1,h2,h3,h4) [diff(h1)-(-max_eig_real*h1+h4);
                   diff(h2)-(-max_eig_real*h2+h3-2*h4);
                   diff(h3)-(-max_eig_real*h3-h1+(2*nu*phi-3*phi^2-mu)*h1);
                   diff(h4)-(-max_eig_real*h4+h2)];
N.lbc = @(h1,h2,h3,h4) [h1+v_u(1); h2+v_u(2);h3+v_u(3);h4+v_u(4)];

[h1,h2,h3,h4] = N\0;

figure(1)
try
    plot3(h1,h2,h3), view(-5,9), %axis off
catch
    plot3(real(h1),imag(h1),abs(h2)), view(-5,9), %axis off
end

n = length(chebcoeffs(h1));

h1_coeffs=chebcoeffs(h1);
h2_coeffs=chebcoeffs(h2);
h3_coeffs=chebcoeffs(h3);
h4_coeffs=chebcoeffs(h4);

psi_coeffs=chebcoeffs(phi);

figure(2)
plot(log(abs(chebcoeffs(h1)))/log(10),'o')
title('Log-Norm of Chebyshev coeff')