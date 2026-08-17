function b1x = get_b1x(params,bndl,mflds,sig,pulse4D,tildes,S)

%sig = get_sig_afterBVP(params,pulse4D,x); %get sigma(x)

[V1s,V2s,V1u,V2u] = get_Vs(params,bndl,mflds,pulse4D,sig(1),sig(2),x,S); %get each soln V(x)

%compute U_1 and U_{\varphi'} using eqn 4.6 from paper3

tbeta = tildes(:,1);
tgamma = tildes(:,2);

bx = tbeta(1)*V1s + tbeta(2)*V2s + tgamma(1)*V1u + tgamma(2)*V2u;

bx = real(bx);

b1x = bx(1);

end