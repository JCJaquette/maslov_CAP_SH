function a1x = get_a1x(params,bndl,mflds,sig,pulse4D,tildes,x,S)

[V1s,V2s,~,~] = get_Vs(params,bndl,mflds,pulse4D,sig(1),sig(2),x,S); %get each soln V(x)

teta = tildes(:,3);

ax = teta(1)*V1s + teta(2)*V2s;

a1x = ax(1);

end