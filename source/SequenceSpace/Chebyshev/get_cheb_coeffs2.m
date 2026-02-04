function seqOut = get_cheb_coeffs2(f_vals,ord)
% take in equispaced function values as row vec, put out cheb coeffs column

    k = chebfun(f_vals','equi');
    seqOut = k.funs{1,1}.onefun.coeffs;
    seqOut(2:end) = seqOut(2:end)/2;
    seqOut = [seqOut;zeros(ord-length(seqOut),1)];

end

