function [V1s,V2s,V1u,V2u] = get_Vs(params,bndl,mflds,sig1,sig2,x,S)

%Start by finding \tilde{V}(x), see eqns below lemma 5.8

tV1s = [exp(mflds.values.s(1) *x);0;0;0];
tV2s = [0;exp(mflds.values.s(2) *x);0;0];

a_res = [bndl.normalForm(1,3,3,1), bndl.normalForm(1,4,2,2);
         bndl.normalForm(2,3,2,2), bndl.normalForm(2,4,1,3)];

tV1u = [a_res(1,1)*x * sig1^2 * exp(mflds.values.s(1) *x);
        a_res(2,1)*x * sig1*sig2 * exp(mflds.values.s(2) *x);
        exp(mflds.values.u(1) *x);
        0];
tV2u = [a_res(1,2)*x * sig1*sig2 * exp(mflds.values.s(1) *x);
        a_res(2,2)*x * sig2^2 * exp(mflds.values.s(2) *x);
        0;
        exp(mflds.values.u(2) *x)];

W_sig = bndl_one_point(real(sig1),imag(sig1),bndl,params);
W_sig = W_sig + infsup(-bndl.r_min.sup,bndl.r_min.sup);

%get V(x) from \tilde{V}(x)

V1s = S * (W_sig*tV1s);
V2s = S * (W_sig*tV2s);
V1u = S * (W_sig*tV1u);
V2u = S * (W_sig*tV2u);

end