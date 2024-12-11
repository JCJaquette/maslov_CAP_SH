clear 

load('verifiedpulse1.mat')
params = getparamsInt(1);

phi_cheb = new_y.a1;

yo1 = (chebcoeff_to_function(new_y.a1))';
yo2 = (chebcoeff_to_function(new_y.a2))';
yo3 = (chebcoeff_to_function(new_y.a3))';
yo4 = (chebcoeff_to_function(new_y.a4))';
phi = [yo1;yo2;yo3;yo4];

phiPrime = RHSofIVP(phi,yo1,params.mu,params.nu,params.L);

%%

load('varbs.mat')

phi = chebfun(1);
phi.domain = [-1,1];
phi.funs{1,1}.onefun.coeffs = [phi_cheb(1),2*phi_cheb(2:end)]';



for i = 1:2

% chebfun integrator

if i == 1
    v_u = real(evecs.u(:,1));
else
    v_u = imag(evecs.u(:,1));
end


N = chebop(-1,1);
N.op = @(t,h1,h2,h3,h4) [diff(h1)-params.L*(h4);
                   diff(h2)-params.L*(h3-2*h4);
                   diff(h3)-params.L*(-h1+(2*params.nu*phi-3*phi^2-params.mu)*h1);
                   diff(h4)-params.L*(h2)];

N.lbc = v_u;
[h1,h2,h3,h4] = N\0;


% newton

n = length(h1);
NN = n;
phi_cheb = [phi_cheb, zeros(1,NN)];
phi_cheb = phi_cheb(1:NN);

h = zeros(1,4*NN);
h(1:n) = chebcoeffs(h1)/2;
h(1) = h(1)*2;
h(NN+1:NN+n) = chebcoeffs(h2)/2;
h(NN+1) = h(NN+1)*2;
h(2*NN+1:2*NN+n) = chebcoeffs(h3)/2;
h(2*NN+1) = h(2*NN+1)*2;
h(3*NN+1:3*NN+n) = chebcoeffs(h4)/2;
h(3*NN+1) = h(3*NN+1)*2;


for j = 1:5

    h = h - (chebDF(phi_cheb,NN,params)\chebF(h,v_u,phi_cheb,NN,params))';

end

disp('norm of F(h) at end of Newton:')
disp(norm(chebF(h,v_u,phi_cheb,NN,params)))

if i == 1                                                                                                                                                                        %hi
    H.a11 = chebcoeff_to_function(h(1:NN));
    H.a21 = chebcoeff_to_function(h(NN+1:2*NN));
    H.a31 = chebcoeff_to_function(h(2*NN+1:3*NN));
    H.a41 = chebcoeff_to_function(h(3*NN+1:4*NN));
else
    H.a12 = chebcoeff_to_function(h(1:NN));
    H.a22 = chebcoeff_to_function(h(NN+1:2*NN));
    H.a32 = chebcoeff_to_function(h(2*NN+1:3*NN));
    H.a42 = chebcoeff_to_function(h(3*NN+1:4*NN));
end


end



detA = (H.a11 .* H.a42) - (H.a12 .* H.a41);
%plot(params.L*(-1:.05:1),detA)





