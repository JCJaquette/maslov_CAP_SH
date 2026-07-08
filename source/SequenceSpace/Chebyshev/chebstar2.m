% the coefficients of order N are stored in the N+1 entry
%
% Chebyshev (cosine) convolution.  Mathematically
%
%     vec(k+1) = sum_{i=-(order-1)}^{order-1}  a(|i|+1) * b(|k-i|+1)
%
% for k = 0..order-1, where terms with |k-i| > order-1 are dropped.
%
% Vectorized over the inner index i.  The previous implementation used a
% scalar double loop (~2*order^2 operations) and a growing array; with INTLAB
% intval arguments that made this routine the single dominant cost of the
% rigorous homoclinic-orbit validation (it is called with order = 3m).
% Vectorizing replaces ~2*order^2 scalar INTLAB operations by ~4*order
% vector ones.  Note sum() uses pairwise summation, so the double-precision
% result can differ from the old sequential sum in the last bits; the intval
% result is rigorous regardless of summation order.
function vec = chebstar2(a,b,order)
    N=max(size(a));
    M=max(size(b));

    if N < order
        a=[a,zeros(1,order-N)];
    end
    if M < order
        b=[b,zeros(1,order-M)];
    end

    i     = -(order-1):(order-1);   % inner summation index
    absi1 = abs(i)+1;
    ai    = a(absi1);               % gather a(|i|+1) once

    vec = (a(1)*b(1))*zeros(1,order);   % zeros of the right type (double/intval)

    for k = 0:order-1
        l   = k - i;
        msk = abs(l) <= order-1;
        vec(k+1) = sum( ai(msk) .* b(abs(l(msk))+1) );
    end

end
