% GETFUNCTIONFROMFOURIERCOEFFS  Computes function values from Fourier
% coefficients on either [-L,L] or [0,L]
%   sol = getFunctionFromFourierCoeffs(S,coeffs, interval_type)
%   sol = S.getFunctionFromFourierCoeffs(coeffs, interval_type)
% 
% returns a matrix of size 2 x N where the first column is the vector of
% times and the second column is the vector of function values. 
function sol = getFunctionFromFourierCoeffs(S,coeffs, interval_type)
order = S.fourier.order; 

if interval_type == "half"
    T=0:.05:S.time;
else
    T=-S.time:.05:S.time;
end

f       = 0*T;
f_d     = 0*T;
f_dd    = 0*T;
f_ddd   = 0*T;
for n = -order:order
    f       = f     +coeffs(n+order+1)*exp(1i*n*pi.*T./S.time);
    f_d     = f_d   +coeffs(n+order+1)*exp(1i*n*pi.*T./S.time)*(1i*n*pi/S.time)^1;
    f_dd    = f_dd  +coeffs(n+order+1)*exp(1i*n*pi.*T./S.time)*(1i*n*pi/S.time)^2;
    f_ddd   = f_ddd +coeffs(n+order+1)*exp(1i*n*pi.*T./S.time)*(1i*n*pi/S.time)^3;
end
f=real(f);
f_d=real(f_d);
f_dd=real(f_dd);
f_ddd=real(f_ddd);

sol(:,1) = T';
sol(:,2) = f';
sol(:,3) = f_d';
sol(:,4) = f_dd';
sol(:,5) = f_ddd';

end
