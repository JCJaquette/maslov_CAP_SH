close all
clear

[params, phi, lambda, mani_coeffs] = getparamsBefore(1);

global_min = exp(-params.new_L);%working with s = exp(t)
global_max = 1;
%domINT = infsup(global_min,global_max);

[dcoeffs1, dcoeffs2] = differentiate_mani(mani_coeffs);
mani_order = length(mani_coeffs) - 1;

sig1 = @(s) s^lambda(1) * phi(1);
sig2 = @(s) s^lambda(2) * phi(2);

detA_taylor = cauchyProd2D(dcoeffs1(:,:,1), dcoeffs2(:,:,4)) ...
                - cauchyProd2D(dcoeffs1(:,:,4), dcoeffs2(:,:,1));
detAPrime_taylor = cauchyProd2D(dcoeffs1(:,:,1), dcoeffs2(:,:,2)) ...
                - cauchyProd2D(dcoeffs1(:,:,2), dcoeffs2(:,:,1));

f = @(s) taylorSum2D(detA_taylor,sig1(s),sig2(s));
df = @(s) taylorSum2D(detAPrime_taylor,sig1(s),sig2(s));

z = linspace(global_min,global_max,500);
plot(z,f(z),'Color','black')
hold on
count0s = 0;
flag = intval(1)*[];
tol = 1e-5;

while isempty(domINT) == 0

    b = domINT(end);

    if f(b.inf)*f(b.sup) > 0

        if abs(f(b) + f_error) > 0

            domINT(end) = [];
            plot([b.inf,b.sup],[0,0],'color','red')

        else

            domINT(end) = [];
            end1 = infsup(b.inf,b.mid);
            end2 = infsup(b.mid,b.sup);
            domINT = [domINT, end1,end2];

        end

    else
        
        if abs(df(b) + df_error) > 0

            domINT(end) = [];
            plot([b.inf,b.sup],[0,0],'color','green')
            count0s = count0s+1;

        else

            domINT(end) = [];
            end1 = infsup(b.inf,b.mid);
            end2 = infsup(b.mid,b.sup);
            domINT = [domINT, end1,end2];

        end

    end

    if rad(b)<tol
        domINT(end) = [];
        flag = [flag;b];
    end

end

if isempty(flag) == 0
    disp('Not validated, check flag')
    return
end
count0s