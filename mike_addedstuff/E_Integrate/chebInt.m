clear 
close all

load('verifiedpulse1.mat')
[params,IC] = getparamsInt(1);
ord = params.ord;

Q = [1, 0, 0, 0;
     0, 0, 1, 0;
     0, 2, 0, 1;
     0, 1, 0, 0]; %Q brings you from thesis Ch9 coords to Ch12 coords

yo1 = (chebcoeff_to_function(new_y.a1))';
yo2 = (chebcoeff_to_function(new_y.a2))';
yo3 = (chebcoeff_to_function(new_y.a3))';
yo4 = (chebcoeff_to_function(new_y.a4))';

pulse9 = [yo1;yo2;yo3;yo4];
pulse12 = Q*pulse9;

phi = pulse9(1,:);
phi_cheb = new_y.a1;% pulse is same in both coordinate systems
phi_cheb = [phi_cheb , zeros(1,ord-length(phi_cheb))];

pulsePrime9 = RHSofODE(pulse9,params.mu,params.nu);
pulsePrime12 = Q*pulsePrime9;
params.cheb.order = ord;
a = get_cheb_coeffs(pulsePrime12, params);
a(101:end,:) = zeros(ord - 100,4);
pulse12ICvec = pulsePrime12(:,1);

h = [a(:,1);a(:,2);a(:,3);a(:,4)]';

    for j = 1:5

        h = h - (chebDF(phi_cheb,ord,params)\chebF(h,pulse12ICvec,phi_cheb,ord,params))';

    end

domn = linspace(-params.L,params.L,201);
hold on
plot(domn,phi)

hold on
plot(domn,pulsePrime12(1,:))

hold on
pulseprimenew1 = chebcoeff_to_function(h(1:ord));
plot(domn,pulseprimenew1)


%%

Ad_N = chebDF(phi_cheb,ord,params);
A_N = inv(Ad_N);

a_bar = h;



%%

phi = chebfun(1);
phi.domain = [-1,1];
phi.funs{1,1}.onefun.coeffs = [phi_cheb(1),2*phi_cheb(2:end)]';

unstableVec1 = zeros(4,1);
unstableVec2 = unstableVec1;

unstableVec1(:,1) = real(IC(1,1,:));
unstableVec2(:,1) = imag(IC(1,1,:));

intICvec = getICvec(pulse12ICvec,unstableVec1,unstableVec2);

N = chebop(-1,1);
N.op = @(t,h1,h2,h3,h4) [diff(h1)-params.L*(h4);
                   diff(h2)-params.L*(h3-2*h4);
                   diff(h3)-params.L*(-h1+(2*params.nu*phi-3*phi^2-params.mu)*h1);
                   diff(h4)-params.L*(h2)];

N.lbc = intICvec;
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

for j = 1:10

    h = h - (chebDF(phi_cheb,NN,params)\chebF(h,intICvec,phi_cheb,NN,params))';

end

disp('norm of F(h) at end of Newton:')
disp(norm(chebF(h,intICvec,phi_cheb,NN,params)))

plot(h1)
hold on
plot(linspace(-1,1,201),chebcoeff_to_function(h(((1-1)*n)+1:1*n)))


%%


H.a11 = chebcoeff_to_function(h(1:NN));
H.a21 = chebcoeff_to_function(h(NN+1:2*NN));
H.a31 = chebcoeff_to_function(h(2*NN+1:3*NN));
H.a41 = chebcoeff_to_function(h(3*NN+1:4*NN));

H.a12 = chebcoeff_to_function(h(1:NN));
H.a22 = chebcoeff_to_function(h(NN+1:2*NN));
H.a32 = chebcoeff_to_function(h(2*NN+1:3*NN));
H.a42 = chebcoeff_to_function(h(3*NN+1:4*NN));


%detA = (H.a11 .* H.a42) - (H.a12 .* H.a41);
%plot(params.L*(-1:.05:1),detA)


%%

Ad_N = chebDF(phi_cheb,NN,params);
A_N = Ad_N^-1;
a_bar = h;
delta = 1.001;

% Y0 = computeY0(A_N,a_bar,phi_cheb,params,NN,v_u,delta)
% 
% Z0 = computeZ0(A_N,Ad_N,NN,delta)

% Ad=intval(Ad);
% A=intval(A);
% a_bar=intval(a_bar);
% delta=intval(delta);
% 
Y0 = computeY0(A,a_bar,phi_cheb,params,NN,ICvec,delta)
% 
Z0 = computeZ0(A,Ad,NN,delta)


