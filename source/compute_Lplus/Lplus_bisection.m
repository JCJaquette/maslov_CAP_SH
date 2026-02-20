function Lplus = Lplus_bisection(fn,eps0V,sigmin)
%given functions fn, and eps0 find Lplus such that eps0V<1 and fn<sigmin

L_int = infsup(1,2);

while fn(L_int.sup) >= sigmin || eps0V(L_int.sup) >= 1
    L_int = infsup(L_int.sup,2*L_int.sup);
end

while L_int.rad > .5
    if fn(L_int.mid) >= sigmin || eps0V(L_int.mid) >= 1
        L_int = infsup(L_int.mid,L_int.sup);
    else
        L_int = infsup(L_int.inf,L_int.mid);
    end
end

Lplus = L_int.sup;

%When i tested this with fn = @(x) 10/x, eps0V = @(x) 20*exp(-x), 
% sigmin = .5, this returns 20, which seems good.

end