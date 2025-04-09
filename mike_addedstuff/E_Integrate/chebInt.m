clear 
close all

load('verifiedpulse2.mat')
[params,mfld_u] = getparamsInt(2);
ord = params.cheb.order;
rho = params.rho;

Q = [1, 0, 0, 0;
     0, 0, 1, 0;
     0, 2, 0, 1;
     0, 1, 0, 0]; %Q brings you from natural coords to skew symmetric coords

x1 = (chebcoeff_to_function(new_y.a1))';
x2 = (chebcoeff_to_function(new_y.a2))';
x3 = (chebcoeff_to_function(new_y.a3))';
x4 = (chebcoeff_to_function(new_y.a4))';

pulse_natural = [x1;x2;x3;x4];
pulse_skewSym = Q*pulse_natural;

phi = pulse_natural(1,:);
phi_cheb = new_y.a1;
phi_cheb = [phi_cheb , zeros(1,ord-length(phi_cheb))];
phi_cheb = phi_cheb(1:ord);

pulsePrime_natural = RHSofODE(pulse_natural,params.mu,params.nu);
pulsePrime_skewSym = Q*pulsePrime_natural;
pulsePrime_skewSym_cheb = get_cheb_coeffs(pulsePrime_skewSym, params);
pulsePrime_skewSym_cheb(101:end,:) = zeros(ord - 100,4);
pulsePrime_skewSym_ICvec = pulsePrime_skewSym(:,1);

h = [pulsePrime_skewSym_cheb(:,1);
    pulsePrime_skewSym_cheb(:,2);
    pulsePrime_skewSym_cheb(:,3);
    pulsePrime_skewSym_cheb(:,4)]';

    for j = 1:5

        h = h - (chebDF(phi_cheb,ord,params)\chebF(h,pulsePrime_skewSym_ICvec,phi_cheb,ord,params))';

    end

% domn = linspace(-params.L,params.L,201);
% hold on
% plot(domn,phi)




%%

phi = chebfun(1);
phi.domain = [-1,1];
phi.funs{1,1}.onefun.coeffs = [phi_cheb(1),2*phi_cheb(2:end)]';

[unstableVec1,unstableVec2] = getEu_minusL(mfld_u.coeffs,mfld_u.pulseIC_phi);
unstableVec_re = real(unstableVec1);
unstableVec_im = imag(unstableVec1);

% transform to symplectic coord 

unstableVec_re_sym = Q*unstableVec_re;
unstableVec_im_sym = Q*unstableVec_im;

intICvec = getICvec(pulsePrime_skewSym_ICvec,unstableVec_re_sym,unstableVec_im_sym);

% plot_manifold(mfld_u.coeffs,25,'red');
% hold on
% k1 = [pulse_natural(1:2,1), pulse_natural(1:2,1)+unstableVec_re(1:2);
%       pulse_natural(4,1), pulse_natural(4,1)+unstableVec_re(4)]';
% k2 = [pulse_natural(1:2,1), pulse_natural(1:2,1)+unstableVec_im(1:2);
%       pulse_natural(4,1), pulse_natural(4,1)+unstableVec_im(4)]';
% pulseprimeicNAT = [pulse_natural(1:2,1), pulse_natural(1:2,1)-pulsePrime_natural(1:2,1);
%       pulse_natural(4,1), pulse_natural(4,1)-pulsePrime_natural(4,1)]';
% intICvecNAT = Q\intICvec;
% intIC = [pulse_natural(1:2,1), pulse_natural(1:2,1)-intICvecNAT(1:2,1);
%       pulse_natural(4,1), pulse_natural(4,1)-intICvecNAT(4,1)]';
% plot3(k1(:,1),k1(:,2),k1(:,3),'black')
% plot3(k2(:,1),k2(:,2),k2(:,3),'black')
% plot3(pulseprimeicNAT(:,1),pulseprimeicNAT(:,2),pulseprimeicNAT(:,3),'black')
% plot3(intIC(:,1),intIC(:,2),intIC(:,3),'black')

Ch12ODE = chebop(-1,1);
Ch12ODE.op = @(t,h1,h2,h3,h4) [diff(h1)-params.L*(h4);
                   diff(h2)-params.L*(h3-2*h4);
                   diff(h3)-params.L*(-h1+(2*params.nu*phi-3*phi^2-params.mu)*h1);
                   diff(h4)-params.L*(h2)];

Ch12ODE.lbc = intICvec;
[h1,h2,h3,h4] = Ch12ODE\0; %#ok<RHSFN>



% newton

n = length(h1);
ord = 600;
phi_cheb = [phi_cheb, zeros(1,ord)];
phi_cheb = phi_cheb(1:ord);

h = zeros(1,4*ord);
h(1:n) = chebcoeffs(h1)/2;
h(1) = h(1)*2;
h(ord+1:ord+n) = chebcoeffs(h2)/2;
h(ord+1) = h(ord+1)*2;
h(2*ord+1:2*ord+n) = chebcoeffs(h3)/2;
h(2*ord+1) = h(2*ord+1)*2;
h(3*ord+1:3*ord+n) = chebcoeffs(h4)/2;
h(3*ord+1) = h(3*ord+1)*2;

for j = 1:3

    h = h - (chebDF(phi_cheb,ord,params)\chebF(h,intICvec,phi_cheb,ord,params))';

end

disp('norm of F(h) at end of Newton:')
disp(norm(chebF(h,intICvec,phi_cheb,ord,params)))

% i = 1; plot(h1)
% hold on
% plot(linspace(-1,1,201),chebcoeff_to_function(h(((i-1)*ord)+1:i*ord)))

% plot(log(abs(h((i-1)*ord+1:i*ord))))


%% CAP

phi_cheb_int = intval(1)*phi_cheb;
Ad_N = chebDF_intval(phi_cheb_int,ord,params);
A_N = Ad_N^-1;
a_bar = intval(1)*h;

%%

Y0 = computeY0(A_N,a_bar,phi_cheb_int,params,ord,intICvec,params.del);%Y0 for pulse 3 is 7.454173333836001e-06
%Y0 = 7.454173333836001e-06;
Y0hat = computeY0hat(A_N,rho,params.L,params.del,params.nu,a_bar(1:600),phi_cheb);
Z0 = computeZ0(A_N,Ad_N,ord,params.del);
Z1 = computeZ1(A_N,ord,phi_cheb_int,params.del,params);
Z2hat = computeZ2hat(A_N,params.L,params.del,params.nu,rho);

rs = 0:10^-12:10^-6;
radii_poly = Y0 + Y0hat - (1-Z0-Z1-Z2hat)*rs;

good_r = sup((Y0 + Y0hat)/(1-Z0-Z1-Z2hat))

% pulse1 soln validated with good_r = 2.4672e-07
% pulse2 soln validated with good_r = 2.4366e-07
% pulse3 soln validated by leaving chebstar2 without fft in y3tail of Y0 function,
% good_r = 7.8659e-06


%%

detA = h1*chebfun(pulsePrime_skewSym(4,:)','equi') - h4*chebfun(pulsePrime_skewSym(1,:)','equi');

plot(detA)

