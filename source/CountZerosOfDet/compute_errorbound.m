function [ferror, dferror] = compute_errorbound(W1coeffs, W2coeffs, W1Primecoeffs, W2Primecoeffs, params, phi)

nu = max(abs(phi));
W_error = 2*pi/log(1/nu) * params.mfld_error;

ferror_coeffs = abs(W1coeffs(:,:,1)) + abs(W2coeffs(:,:,2)) - ... %These show up
                    abs(W1coeffs(:,:,2)) - abs(W2coeffs(:,:,1));  %in the f error 
ferror_coeffsum = sum(sum(ferror_coeffs));
ferror = W_error*(ferror_coeffsum + 2*eps);

W1prime_error = (1 + abs(params.lambda(1)))*W_error;
W2prime_error = (1 + abs(params.lambda(2)))*W_error;

df_coeffsum11 = sum(sum(abs(W1coeffs(:,:,1)) + abs(W1coeffs(:,:,2))));
df_coeffsum12 = sum(sum(abs(W2coeffs(:,:,1)) + abs(W2coeffs(:,:,2))));
df_coeffs2 = abs(W1Primecoeffs(:,:,1)) + abs(W1Primecoeffs(:,:,2)) + ...
                 abs(W2Primecoeffs(:,:,1)) + abs(W2Primecoeffs(:,:,2));
df_coeffsum2 = sum(sum(df_coeffs2));

dferror = W2prime_error*df_coeffsum11 + W1prime_error*df_coeffsum12 + ...
            W_error*df_coeffsum2 + 2*W_error*(W1prime_error + W2prime_error);

end