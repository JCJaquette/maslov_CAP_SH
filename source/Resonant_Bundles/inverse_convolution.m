function [W_coeff_inv] = inverse_convolution(W_coeff,N_inv)
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here

disp('Computing Inverse')

[N,M]=size(W_coeff);
if N~=M
    disp("WARNING NON SQUARE DIMENSION!")
end
N_inv = max([N_inv,N]) ;

N_inv=N_inv+1;

% N_inv=N;
W_coeff_inv = zeros(N_inv,N_inv); 
base =W_coeff(1,1);
try
    W_coeff_inv(1,1)=1/base;
catch
    W_coeff_inv(1,1)=1/mid(base);
end
% W_coeff_inv = ifft(fft(coeff_pad).^-1);
% W_coeff_inv = W_coeff_inv(1:N+1,1:N+1);


W_null = W_coeff;
W_null(1,1)=0;
tic

for alphaNorm=1:(N_inv-1)
    % alphaNorm
    % prod = Cauchyfft2(W_coeff_inv,W_null);
    prod = conv2(W_coeff_inv,W_null);
    for i=0:alphaNorm
        W_coeff_inv(1+i,alphaNorm-i+1) = -prod(1+i,alphaNorm-i+1)/base;
    end
end
toc

prod = conv2(W_coeff_inv,W_coeff);


end