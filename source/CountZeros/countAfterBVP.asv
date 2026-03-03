function [count0s,flag] = countAfterBVP(params,mflds,pulse4D,tol,BOOL_plot)
% Count zeros of the determinant on [-L_conj, -L_bvp]

global_min = 0;
global_max = sup(mflds.Lplus - pulse4D.Lbvp);
domINT = infsup(global_min,global_max);

[W1,W2] = differentiate_mani(mid(mflds.stable.coeffs));
W1 = 1i*W1;

Lambda = get_Lambda(params,'s');

W1Prime = W1(:,:,2:4) - Lambda(1)*W1(:,:,1:3); %See Lemma 5.6
W2Prime = W2(:,:,2:4) - Lambda(2)*W2(:,:,1:3); 

theta = [params.rho*exp(pulse4D.psi*1i), params.rho*exp(pulse4D.psi*1i)];

[f_error, df_error] = get_detW_errorbound(W1, W2, W1Prime, W2Prime, mflds, theta, Lambda);

sig1 = @(t) exp(t*Lambda(1)) * theta(1); %sig(t) = [sig1;sig2] 
sig2 = @(t) exp(t*Lambda(2)) * theta(2); %       = exp(Lambda*t)*theta

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