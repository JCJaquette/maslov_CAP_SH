function Eu = chebInt(params,mflds,y)
% Get initial condition and integrate it

params = struct_intvaltodouble(params);%Standard numerics, we use doubles
mflds = struct_intvaltodouble(mflds);
y = struct_intvaltodouble(y);

mfld_u.coeffs = mflds.unstable.coeffs;
mfld_u.pulseIC_phi = [y.phi1,y.phi2];
ord = params.Eu.order;

Q = [1, 0, 0, 0; 
     0, 0, 1, 0;
     0, 2, 0, 1;
     0, 1, 0, 0]; %Q brings you from natural coords to skew symmetric coords

pulse_natural_cheb = [y.a1; y.a2; y.a3; y.a4];
pulse_natural_cheb = [pulse_natural_cheb, zeros(4,ord - length(pulse_natural_cheb))];
pulse_natural_cheb = pulse_natural_cheb(:,1:ord);
phi_cheb = pulse_natural_cheb(1,:);

pulsePrime_natural_cheb = RHSofODE_coeffs(pulse_natural_cheb,params,y.Lbvp);
Eu.U_vpp_cheb = (Q*pulsePrime_natural_cheb)';


for i = 1:4
    U_vp_ICvec(i) = [chebSum(Eu.U_vpp_cheb(:,i),-1)];
end


phi = chebfun(1);
phi.domain = [-1,1];
phi.funs{1,1}.onefun.coeffs = [phi_cheb(1),2*phi_cheb(2:end)]';

[unstableVec1,~] = getEu_minusL(mfld_u.coeffs,mfld_u.pulseIC_phi);
unstableVec_re = real(unstableVec1);
unstableVec_im = imag(unstableVec1);

% transform to symplectic coordinates

unstableVec_re_sym = Q*unstableVec_re;
unstableVec_im_sym = Q*unstableVec_im;

intICvec = getICvec(U_vp_ICvec,unstableVec_re_sym,unstableVec_im_sym);

ODE = chebop(-1,1);
ODE.op = @(t,h1,h2,h3,h4) [diff(h1)-y.Lbvp*(h4);
                   diff(h2)-y.Lbvp*(h3-2*h4);
                   diff(h3)-y.Lbvp*(-h1+(2*params.nu*phi-3*phi^2-params.mu)*h1);
                   diff(h4)-y.Lbvp*(h2)];

ODE.lbc = intICvec;
[h1,h2,h3,h4] = ODE\0; %#ok<RHSFN>

% Set the coefficients into the form we want(a single 4 x ord matrix)

length_vec = [length(h1);length(h2);length(h3);length(h4)];

n = max(length_vec);
Eu.nonzero = 2^(ceil(log2(n)));

h_cheb = zeros(Eu.nonzero,4);
h_cheb(1:length_vec(1),1) = chebcoeffs(h1)/2; 
h_cheb(1:length_vec(2),2) = chebcoeffs(h2)/2;
h_cheb(1:length_vec(3),3) = chebcoeffs(h3)/2;
h_cheb(1:length_vec(4),4) = chebcoeffs(h4)/2;
h_cheb(1,:) = 2*h_cheb(1,:); 
Eu.U_1_cheb = h_cheb;

end