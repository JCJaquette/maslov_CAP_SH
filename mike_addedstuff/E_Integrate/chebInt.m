clear 
close all

load('verifiedpulse2.mat')
[params,mfld_u] = getparamsInt(2);
mfld_u.pulseIC_phi = [new_y.phi1,new_y.phi2];
ord = params.cheb.order;
rho = params.rho;

Q = [1, 0, 0, 0; 
     0, 0, 1, 0;
     0, 2, 0, 1;
     0, 1, 0, 0]; %Q brings you from natural coords to skew symmetric coords

pulse_natural_cheb = [new_y.a1; new_y.a2; new_y.a3; new_y.a4];
pulse_natural_cheb = [pulse_natural_cheb, zeros(4,ord - length(pulse_natural_cheb))];
pulse_natural_cheb = pulse_natural_cheb(:,1:ord);
phi_cheb = pulse_natural_cheb(1,:);
pulse_skewSym_cheb = Q*pulse_natural_cheb;

pulsePrime_natural_cheb = RHSofODE_coeffs(pulse_natural_cheb,params);
pulsePrime_skewSym_cheb = (Q*pulsePrime_natural_cheb)';


for i = 1:4
    pulsePrime_skewSym_ICvec(i) = [chebSum(pulsePrime_skewSym_cheb(:,i),-1)];
end

h = [pulsePrime_skewSym_cheb(:,1);
    pulsePrime_skewSym_cheb(:,2);
    pulsePrime_skewSym_cheb(:,3);
    pulsePrime_skewSym_cheb(:,4)]';

norm(chebF(h,pulsePrime_skewSym_ICvec,phi_cheb,ord,params),1)


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
% pulse_natural = chebcoeff_to_function(pulse_natural_cheb);
% pulsePrime_natural = chebcoeff_to_function(pulsePrime_natural_cheb);
% % k1 = [pulse_natural(1,1:2), pulse_natural(1,4);
% %       pulse_natural(1,1:2)+unstableVec_re(1:2)', pulse_natural(1,4)+unstableVec_re(4)];
% % k2 = [pulse_natural(1,1:2), pulse_natural(1,4);
% %       pulse_natural(1,1:2)+unstableVec_im(1:2)', pulse_natural(1,4)+unstableVec_im(4)];
% pulseprimeicNAT = [pulse_natural(1,1:2), pulse_natural(1,4);
%       pulse_natural(1,1:2)+pulsePrime_natural(1,1:2)/50, pulse_natural(1,4)+pulsePrime_natural(1,4)/50];
% intICvecNAT = Q\intICvec;
% intIC = [pulse_natural(1,1:2), pulse_natural(1,4);
%       pulse_natural(1,1:2)+intICvecNAT(1:2,1)', pulse_natural(1,4)+intICvecNAT(4,1)];
% % plot3(k1(:,1),k1(:,2),k1(:,3),'black')
% % plot3(k2(:,1),k2(:,2),k2(:,3),'black')
% plot3(pulseprimeicNAT(:,1),pulseprimeicNAT(:,2),pulseprimeicNAT(:,3),'black')
% plot3(intIC(:,1),intIC(:,2),intIC(:,3),'black')
% plot3(pulse_natural(:,1),pulse_natural(:,2),pulse_natural(:,4),'color',[0.4940 0.1840 0.5560])

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

% for j = 1:1
% 
%     h = h - (chebDF(phi_cheb,ord,params)\chebF(h,intICvec,phi_cheb,ord,params))';
% 
% end
% 
% disp('norm of F(h) at end of Newton:')
% disp(norm(chebF(h,intICvec,phi_cheb,ord,params)))

% i = 1; plot(h1)
% hold on
% plot(linspace(-1,1,201),chebcoeff_to_function(h(((i-1)*ord)+1:i*ord)))

% plot(log(abs(h((i-1)*ord+1:i*ord))))

figure;

subplot(4, 1, 1);
plot(h1);
subplot(4, 1, 2);
plot(h2);
subplot(4, 1, 3);
plot(h3);
subplot(4, 1, 4);
plot(h4);

sgtitle('Numerical Solution');

1
%% CAP

phi_cheb_int = intval(1)*phi_cheb;
Ad_N = chebDF_intval(phi_cheb_int,ord,params);
A_N = Ad_N^-1;
a_bar = intval(1)*h;


Y0 = computeY0(A_N,a_bar,phi_cheb_int,params,ord,intICvec,params.del);
Y0hat = computeY0hat(A_N,rho,params.L,params.del,params.nu,a_bar(1:600),phi_cheb);
Z0 = computeZ0(A_N,Ad_N,ord,params.del);
Z1 = computeZ1(A_N,ord,phi_cheb_int,params.del,params);
Z2hat = computeZ2hat(A_N,params.L,params.del,params.nu,rho);

rs = 0:10^-12:10^-6;
radii_poly = Y0 + Y0hat - (1-Z0-Z1-Z2hat)*rs;

good_r = sup((Y0 + Y0hat)/(1-Z0-Z1-Z2hat))

% pulse3 soln validated by leaving chebstar2 without fft in y3tail of Y0 function


%%

%detA = h1*chebfun(pulsePrime_skewSym(4,:)','equi') - h4*chebfun(pulsePrime_skewSym(1,:)','equi');
detA = h1*chebfun(chebcoeff_to_function(pulsePrime_skewSym_cheb(:,4)'),'equi') ...
    - h4*chebfun(chebcoeff_to_function(pulsePrime_skewSym_cheb(:,1)'),'equi');

plot(detA)

%%
m = length(h)/4;
h_cheb = [h(1:m); h(m+1:2*m); h(2*m+1:3*m); h(3*m+1:4*m)]';
phiPrime_cheb = pulsePrime_skewSym_cheb;
phi_cheb = pulse_skewSym_cheb';
%E_h = good_r;
E_h = 0;

save('varbs3','phi_cheb','h_cheb','phiPrime_cheb','E_h')