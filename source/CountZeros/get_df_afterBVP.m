function dfx = get_df_afterBVP(params,bndl,mflds,pulse4D,tildes,x,S)

sig = get_sig_afterBVP(params,pulse4D,x); %get sigma(x)

[V1s,V2s,V1u,V2u] = get_Vs(params,bndl,mflds,sig(1),sig(2),x,S); %get each soln V(x)

%compute U_1 and U_{\varphi'} using eqn 4.6 from paper3

tbeta = tildes(:,1);
tgamma = tildes(:,2);
teta = tildes(:,3);

ax = teta(1)*V1s + teta(2)*V2s;

bx = tbeta(1)*V1s + tbeta(2)*V2s + tgamma(1)*V1u + tgamma(2)*V2u;

dfx = ax(1)*bx(2) - ax(2)*bx(1);

end