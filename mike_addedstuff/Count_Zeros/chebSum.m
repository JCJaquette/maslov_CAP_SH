function out = chebSum(a,x)
%cheb coeffs in, function at x out

Ns = 0:length(a)-1;
a(2:end) = 2*a(2:end);

if size(Ns) ~= size(a)
    Ns = Ns';
end

out = sum(a.*cos(Ns*acos(x)));

end

