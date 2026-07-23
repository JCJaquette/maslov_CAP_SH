function [count0s,flag] = countAfterBVP(params,bndl,mflds,pulse4D,Eu,tol,BOOL_plot)
% Count zeros of the determinant on [-L_conj, -L_bvp]

global_min = 0;
global_max = sup(mflds.Lplus - pulse4D.Lbvp);
domINT = infsup(global_min,global_max);

%get eta,beta,gamma

[tbeta,tgamma] = getLinSolnCoords(params,bndl,Eu.U_1_cheb,pulse4D);

teta = getLinSolnCoords(params,bndl,Eu.U_vpp_cheb,pulse4D);

tildes = [tbeta,tgamma,teta];

S = [1, 0, 0, 0; 
     0, 0, 1, 0;
     0, 2, 0, 1;
     0, 1, 0, 0];

f = @(x) get_F_afterBVP(params,bndl,mflds,pulse4D,tildes,x,S);
df = @(x) get_df_afterBVP(params,bndl,mflds,pulse4D,tildes,x,S);

[count0s, flag] = getZeroCount(domINT, f, 0, df, 0, tol, BOOL_plot);
%note ferror and dferror set to zero since it's all contained in intvals
%^^ TODO: CHECK THIS ^^

if 1

    figure
    testX = linspace(global_min,global_max,1000);
    fX = intval(0)*testX; testdf = fX; numdf = fX;
    dx = testX(2)-testX(1);
    for i = 2:999
        numdf(i) = (f(testX(i+1)) - f(testX(i-1)))/(2*dx);
    end
    for i = 2:999
        testdf(i) = get_df_afterBVP(params,bndl,mflds,pulse4D,tildes,testX(i),S);
    end
    for i = 2:999
        fX(i) = get_F_afterBVP(params,bndl,mflds,pulse4D,tildes,testX(i),S);
    end

    % plot(testX,mid(numdf))
    % hold on
    % plot(testX,mid(testdf))
    % plot(testX,mid(fX))

    plot(testX,numdf)
    hold on
    plot(testX,testdf)
    plot(testX,fX)

%in df function also test for a' and b'
end

end