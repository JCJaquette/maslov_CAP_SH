function prod_trunc = Cauchyfft2( A,B)
N1 = length(A);
N2 = length(B);
N = N1+N2;

coeff_pad = intval(0)*zeros(2*N,2*N);
A_pad = coeff_pad;
B_pad = coeff_pad;
A_pad(1:N1,1:N1)=A;
B_pad(1:N2,1:N2)=B;


% 
prod = (verifyfft(verifyfft(A_pad,1).',1).') .* (verifyfft(verifyfft(B_pad,1).',1).');
prod = verifyfft(verifyfft(prod,-1).',-1).';
prod_trunc = prod(1:N,1:N);

end