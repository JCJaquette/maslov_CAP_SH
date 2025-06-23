close all
clear

[params, phi, lambda, mani_coeffs] = getparamsBefore(1);
M = real(lambda(1));

global_min = -params.new_L;
global_max = 0;
domINT = infsup(global_min,global_max);

[a,b] = differentiate_mani(mani_coeffs);
a = 1i*a;
mani_order = length(mani_coeffs) - 1;

detA_taylor = cauchyProd2D(a(:,:,1), b(:,:,4)) ...
                - cauchyProd2D(a(:,:,4), b(:,:,1));
detAPrime_taylor = cauchyProd2D(a(:,:,1), b(:,:,2)) ...
                - cauchyProd2D(a(:,:,2), b(:,:,1));

sig1 = @(t) exp(t*lambda(1)) * (phi(1) - 1i*phi(2));
sig2 = @(t) exp(t*lambda(2)) * (phi(1) + 1i*phi(2));

f = @(t) taylorSum2D(detA_taylor,sig1(t),sig2(t));
df = @(t) taylorSum2D(detAPrime_taylor,sig1(t),sig2(t)) - 2*M*f(t);

%%
% params.mfld.order = mani_order;
% for k = 1:4
%     ictest(k,1) = taylorSum2D(mani_coeffs(:,:,k),sig1(0),sig2(0));
% end
% mfldic = mfld_one_point(phi(1),phi(2),mani_coeffs,params);
% load('varbs1.mat')
% varphi_ic_check1 = 0;
% varphi_ic_check2 = 0;
% for k = 1:4
%     varphik = chebcoeff_to_function(phi_cheb(:,k)');
%     varphi_ic_check1(k,1) = varphik(1);
%     varphi_ic_check2(k,1) = chebSum(phi_cheb(:,k),-1);
% end

%%
z = linspace(global_min,global_max,500);
for k = 1:500
    fz(k) = real(f(z(k)));
end
plot(z,fz,'Color','blue')
hold on
for k = 1:500
    dfz(k) = real(df(z(k)));
end
% plot(z,dfz,'Color','red')
plot(z,0*z,'color','black')
roots = FindZero1D(z,dfz);
for i = 1:length(roots)
plot([roots(i),roots(i)],1.5*[min([fz,dfz]),max([fz,dfz])],'color','green')
end


%%
count0s = 0;
flag = intval(1)*[];
tol = 1e-5;
f_error = 0; df_error = 0;

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