function Eu = chebInt(params,mflds,y,Euplot)
% Get initial condition and integrate it

params = struct_intvaltodouble(params);%Standard numerics, we use doubles
mflds = struct_intvaltodouble(mflds);
y = struct_intvaltodouble(y);

mfld_u.coeffs = mflds.unstable.coeffs;
mfld_u.pulseIC_phi = [y.phi1,y.phi2];
ord = params.Eu.order;

S = [1, 0, 0, 0; 
     0, 0, 1, 0;
     0, 2, 0, 1;
     0, 1, 0, 0]; %S brings you from natural coords to symplectic coords

pulse_natural_cheb = [y.a1; y.a2; y.a3; y.a4];
pulse_natural_cheb = [pulse_natural_cheb, zeros(4,ord - length(pulse_natural_cheb))];
pulse_natural_cheb = pulse_natural_cheb(:,1:ord);
phi_cheb = pulse_natural_cheb(1,:);

pulsePrime_natural_cheb = RHSofODE_coeffs(pulse_natural_cheb,params,y.Lbvp);
Eu.U_vpp_cheb = (S*pulsePrime_natural_cheb)';
Eu.U_vpp_r = S*[y.r; y.r; y.r; 
              (params.mu + 3)*y.r + 4*params.nu*y.r^2 + 16*y.r^3];


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

unstableVec_re_sym = S*unstableVec_re;
unstableVec_im_sym = S*unstableVec_im;

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
Eu.U_vpp_r = max(Eu.U_vpp_r' + sum(abs(Eu.U_vpp_cheb(n+1:end,:))));
Eu.U_vpp_cheb = Eu.U_vpp_cheb(1:Eu.nonzero,:);

if Euplot

    labels = {'h_1','h_2','h_3','h_4'};
    figure(1); clf;
    for i = 1:4
        c = Eu.U_vpp_cheb(:,i);
        c = c(:);
    
        f = chebfun(1);
        f.domain = [-1, 1];
        f.funs{1,1}.onefun.coeffs = [c(1); 2*c(2:end)];
        f = newDomain(f, [-y.Lbvp, y.Lbvp]);
    
        subplot(2,2,i);
        plot(f, 'LineWidth', 1.5);
        title(labels{i});
        grid on;
    end
    %sgtitle('\varphi');
    
    % --- plot U_1 components (integrated solution, plotted directly) ---
    hs = {h1, h2, h3, h4};
    figure(2); clf;
    for i = 1:4
        hp = newDomain(hs{i}, [-y.Lbvp, y.Lbvp]);
    
        subplot(2,2,i);
        plot(hp, 'LineWidth', 1.5);
        title(labels{i});
        grid on;
    end
    %sgtitle('U_1  (integrated homogeneous solution)');
    
labels = {'h_1','h_2','h_3','h_4'};
hs = {h1, h2, h3, h4};

figure(3); clf;
for i = 1:4
    % --- left column: chebyshev/phi solution ---
    c = Eu.U_vpp_cheb(:,i);
    c = c(:);
    f = chebfun(1);
    f.domain = [-1, 1];
    f.funs{1,1}.onefun.coeffs = [c(1); 2*c(2:end)];
    f = newDomain(f, [-y.Lbvp, y.Lbvp]);

    subplot(4,2,2*(i-1)+1);
    plot(f, 'LineWidth', 1.5);
    %title([labels{i} ' (\varphi)']);
    grid on;

    % --- right column: integrated homogeneous solution ---
    hp = newDomain(hs{i}, [-y.Lbvp, y.Lbvp]);

    subplot(4,2,2*(i-1)+2);
    plot(hp, 'LineWidth', 1.5);
    %title([labels{i} ' (U_1)']);
    grid on;
end
%sgtitle('\varphi (left) vs U_1 (right)');

end


end