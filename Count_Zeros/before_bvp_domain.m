% Count zeros of the determinant on [-L_conj, -L_bvp]

close all
clear

[params, phi, mani_coeffs] = getparamsBefore(2);

global_min = -params.new_L;
global_max = 0;
domINT = infsup(global_min,global_max);

[W1,W2] = differentiate_mani(mani_coeffs);
W1 = 1i*W1;

W1Prime = W1(:,:,2:4) - params.lambda(1)*W1(:,:,1:3); 
W2Prime = W2(:,:,2:4) - params.lambda(2)*W2(:,:,1:3); 

[f_error, df_error] = compute_errorbound(W1, W2, W1Prime, W2Prime, params, phi);

sig1 = @(t) exp(t*params.lambda(1)) * (phi(1) + phi(2)*1i);
sig2 = @(t) exp(t*params.lambda(2)) * (phi(1) - phi(2)*1i);

detA_taylor = cauchyProd2D(W1(:,:,1), W2(:,:,2)) ...
                - cauchyProd2D(W1(:,:,2), W2(:,:,1));

detAPrime_taylor = cauchyProd2D(W1Prime(:,:,1), W2(:,:,2)) ...
                  + cauchyProd2D(W1(:,:,1), W2Prime(:,:,2)) ...
                  - cauchyProd2D(W1Prime(:,:,2), W2(:,:,1)) ...
                  - cauchyProd2D(W1(:,:,2), W2Prime(:,:,1));

f = @(t) taylorSum2D(detA_taylor,sig1(t),sig2(t));
df = @(t) taylorSum2D(detAPrime_taylor,sig1(t),sig2(t));
numdf = @(t) ( f(t+10^-6) - f(t))/10^-6;

%%

% params.mfld.order = length(mani_coeffs) - 1;
% for k = 1:4
%     ictest(k,1) = taylorSum2D(mani_coeffs(:,:,k),sig1(0),sig2(0));
% end
% mfldic = mfld_one_point(phi(1),phi(2),mani_coeffs,params);
% load('varbs1.mat')
% varphi_ic_check1 = 0;
% varphi_ic_check2 = 0;
% for k = 1:4
%     varphik = chebcoeff_to_function(phi_cheb(:,k)');
%     varphi_ic_check1(k,1) = varphik(1);
%     varphi_ic_check2(k,1) = chebSum(phi_cheb(:,k),-1);
% end


%%

% hold on
% z = linspace(global_min,global_max,500);
% for k = 1:500
%     fz(k) = real(f(z(k)));
% end
% plot(z,fz,'Color','black')
% for k = 1: 500
%     dfz(k) = real(df(z(k)));
% end
% plot(z,dfz,'Color','red')
% for k = 1:500
%     numdfz(k) = real(numdf(z(k)));
% end
% plot(z,numdfz,'color','green')
% roots = FindZero1D(z,dfz);
% for i = 1:length(roots)
% plot([roots(i),roots(i)],1.5*[min([fz,dfz]),max([fz,dfz])],'color','green')
% end

%%

tol = 1e-5;

[count0s, flag] = getZeroCount(domINT, f, f_error, df, df_error, tol);

count0s