function [U_vpp_cheb,U_1_cheb] = chebInt(params,mflds,y)

mfld_u.coeffs = mflds.unstable.coeffs;
mfld_u.pulseIC_phi = [y.phi1,y.phi2];
ord = params.cheb.order;
rho = params.rho;

Q = [1, 0, 0, 0; 
     0, 0, 1, 0;
     0, 2, 0, 1;
     0, 1, 0, 0]; %Q brings you from natural coords to skew symmetric coords

pulse_natural_cheb = [y.a1; y.a2; y.a3; y.a4];
pulse_natural_cheb = [pulse_natural_cheb, zeros(4,ord - length(pulse_natural_cheb))];
pulse_natural_cheb = pulse_natural_cheb(:,1:ord);
phi_cheb = pulse_natural_cheb(1,:);
pulse_skewSym_cheb = Q*pulse_natural_cheb;

pulsePrime_natural_cheb = RHSofODE_coeffs(pulse_natural_cheb,params);
U_vpp_cheb = (Q*pulsePrime_natural_cheb)';


for i = 1:4
    U_vp_ICvec(i) = [chebSum(U_vpp_cheb(:,i),-1)];
end


% h = [U_vp_cheb(:,1);
%     U_vp_cheb(:,2);
%     U_vp_cheb(:,3);
%     U_vp_cheb(:,4)]';
% 
% norm(chebF(h,U_vp_ICvec,phi_cheb,ord,params),1)


%

phi = chebfun(1);
phi.domain = [-1,1];
phi.funs{1,1}.onefun.coeffs = [phi_cheb(1),2*phi_cheb(2:end)]';

[unstableVec1,~] = getEu_minusL(mfld_u.coeffs,mfld_u.pulseIC_phi);
unstableVec_re = real(unstableVec1);
unstableVec_im = imag(unstableVec1);

% transform to symplectic coord 

unstableVec_re_sym = Q*unstableVec_re;
unstableVec_im_sym = Q*unstableVec_im;

intICvec = getICvec(U_vp_ICvec,unstableVec_re_sym,unstableVec_im_sym);

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
Ch12ODE.op = @(t,h1,h2,h3,h4) [diff(h1)-params.Lbvp*(h4);
                   diff(h2)-params.Lbvp*(h3-2*h4);
                   diff(h3)-params.Lbvp*(-h1+(2*params.nu*phi-3*phi^2-params.mu)*h1);
                   diff(h4)-params.Lbvp*(h2)];

Ch12ODE.lbc = intICvec;
[h1,h2,h3,h4] = Ch12ODE\0; %#ok<RHSFN>

% Set the coefficients into the form we want(a single 4 x ord matrix)

n = length(h1);
ord = 600;
phi_cheb = [phi_cheb, zeros(1,ord)];
phi_cheb = phi_cheb(1:ord);

U_1_cheb = zeros(ord,4);
U_1_cheb(1:n,1) = chebcoeffs(h1)/2;
U_1_cheb(1:n,2) = chebcoeffs(h2)/2;
U_1_cheb(1:n,3) = chebcoeffs(h3)/2;
U_1_cheb(1:n,4) = chebcoeffs(h4)/2;
U_1_cheb(1,:) = 2*U_1_cheb(1,:);

% U_1_cheb = zeros(1,4*ord);
% U_1_cheb(1:n) = chebcoeffs(h1)/2;
% U_1_cheb(1) = U_1_cheb(1)*2;
% U_1_cheb(ord+1:ord+n) = chebcoeffs(h2)/2;
% U_1_cheb(ord+1) = U_1_cheb(ord+1)*2;
% U_1_cheb(2*ord+1:2*ord+n) = chebcoeffs(h3)/2;
% U_1_cheb(2*ord+1) = U_1_cheb(2*ord+1)*2;
% U_1_cheb(3*ord+1:3*ord+n) = chebcoeffs(h4)/2;
% U_1_cheb(3*ord+1) = U_1_cheb(3*ord+1)*2;

% for j = 1:1
% 
%     h = h - (chebDF(phi_cheb,ord,params)\chebF(h,intICvec,phi_cheb,ord,params))';
% 
% end
% 
% disp('norm of F(h) after Newton:')
% disp(norm(chebF(h,intICvec,phi_cheb,ord,params)))
%
% i = 1; plot(h1)
% hold on
% plot(linspace(-1,1,201),chebcoeff_to_function(h(((i-1)*ord)+1:i*ord)))
%
% plot(log(abs(h((i-1)*ord+1:i*ord))))
%
% figure;
% 
% subplot(4, 1, 1);
% plot(h1);
% subplot(4, 1, 2);
% plot(h2);
% subplot(4, 1, 3);
% plot(h3);
% subplot(4, 1, 4);
% plot(h4);
% 
% sgtitle('Numerical Solution');

m = length(U_1_cheb)/4;
h_cheb = [U_1_cheb(1:m); U_1_cheb(m+1:2*m); U_1_cheb(2*m+1:3*m); U_1_cheb(3*m+1:4*m)]';
phiPrime_cheb = U_vpp_cheb;
phi_cheb = pulse_skewSym_cheb';
% E_h = good_r;
% E_phi = 0;
% 
% save('varbs3','phi_cheb','h_cheb','phiPrime_cheb','E_h','E_phi')


end