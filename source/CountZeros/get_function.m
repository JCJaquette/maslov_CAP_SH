function fx = get_function(params,bndl,mflds,pulse4D,tildes,x)
%todo: figure out error

% Get the solutions V(x)

Lambda = get_Lambda(params,'s');

theta = [params.rho * exp(pulse4D.psi*1i);
         params.rho * exp(-pulse4D.psi*1i)];

sig1 = exp(x*Lambda(1)) * theta(1); %sig(x) = [sig1;sig2] 
sig2 = exp(x*Lambda(2)) * theta(2); %       = exp(Lambda*x)*theta

a_res = [bndl.normalForm(1,3,3,1), bndl.normalForm(1,4,2,2);
               bndl.normalForm(2,3,2,2), bndl.normalForm(2,4,1,3)];

%Start by finding \tilde{V}(x), see eqns below lemma 5.8

tV1s = [exp(mflds.values.s(1) *x);0;0;0];
tV2s = [0;exp(mflds.values.s(2) *x);0;0];

tV1u = [a_res(1,1)*x * sig1^2 * exp(mflds.values.s(1) *x);
        a_res(2,1)*x * sig1*sig2 * exp(mflds.values.s(2) *x);
        exp(mflds.values.u(1) *x);
        0];
tV2u = [a_res(1,2)*x * sig1*sig2 * exp(mflds.values.s(1) *x);
        a_res(2,2)*x * sig2^2 * exp(mflds.values.s(2) *x);
        0;
        exp(mflds.values.u(2) *x)];

W_sig = bndl_one_point(sig1,sig2,bndl,params);

%get V(x) from \tilde{V}(x)

V1s = W_sig*tV1s;
V2s = W_sig*tV2s;
V1u = W_sig*tV1u;
V2u = W_sig*tV2u;

%compute a and b using eqn 6.1 from paper2

tbeta = tildes(:,1);
tgamma = tildes(:,2);
teta = tildes(:,3);

ax = teta(1)*V1s + teta(2)*V2s;

bx = tbeta(1)*V1s + tbeta(2)*V2s + tgamma(1)*V1u + tgamma(2)*V2u;

fx = ax(1)*bx(4) - ax(4)*bx(1);

end