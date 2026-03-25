function [count0s,flag] = countBeforeBVP(params,mflds,pulse4D,tol,BOOL_plot)
% Count zeros of the determinant on [-L_conj, -L_bvp]

global_min = inf(-mflds.Lminus + pulse4D.Lbvp);
global_max = 0;
domINT = infsup(global_min,global_max);

[W1,W2] = differentiate_mani(mid(mflds.unstable.coeffs));
W1 = 1i*W1;

Lambda = get_Lambda(params,'u');

W1Prime = W1(:,:,2:4) - Lambda(1)*W1(:,:,1:3); %See Lemma 5.6
W2Prime = W2(:,:,2:4) - Lambda(2)*W2(:,:,1:3); 

phi = [pulse4D.phi1, pulse4D.phi2];

[f_error, df_error] = get_detW_errorbound(W1, W2, W1Prime, W2Prime, mflds, phi, Lambda);

sig1 = @(t) exp(t*Lambda(1)) * (phi(1) + phi(2)*1i); %sig(t) = [sig1;sig2] 
sig2 = @(t) exp(t*Lambda(2)) * (phi(1) - phi(2)*1i); %       = exp(Lambda*t)*phi

detA_taylor = cauchyProd2D(W1(:,:,1), W2(:,:,2)) ...
                - cauchyProd2D(W1(:,:,2), W2(:,:,1)); %Get det(A), note we use W_i(:,:,1/2) 
                                                      %since W is in natural coords but this
                                                      %calculation is in symplectic

detAPrime_taylor = cauchyProd2D(W1Prime(:,:,1), W2(:,:,2)) ...
                  + cauchyProd2D(W1(:,:,1), W2Prime(:,:,2)) ...
                  - cauchyProd2D(W1Prime(:,:,2), W2(:,:,1)) ...
                  - cauchyProd2D(W1(:,:,2), W2Prime(:,:,1)); %Get its derivative

f = @(t) taylorSum2D(detA_taylor,sig1(t),sig2(t));
df = @(t) taylorSum2D(detAPrime_taylor,sig1(t),sig2(t));

[count0s, flag] = getZeroCount(domINT, f, f_error, df, df_error, tol, BOOL_plot);

end