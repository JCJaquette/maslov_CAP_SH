close all
clear

global_min = -1;
global_max = 1;
domINT = infsup(global_min,global_max);

[params,h_cheb,phi_cheb,phiPrime_cheb] = getparamsCount(3);

f_error = infsup(-params.f_error,params.f_error);
df_error = infsup(-params.df_error,params.df_error);

A_cheb = chebstar2fft(h_cheb(:,1),phiPrime_cheb(:,4)) ...
    - chebstar2fft(h_cheb(:,4),phiPrime_cheb(:,1));
APrime_cheb = chebstar2fft(h_cheb(:,1),phiPrime_cheb(:,2)) ...
    - chebstar2fft(h_cheb(:,2),phiPrime_cheb(:,1));

f = @(x) chebSum(A_cheb,x);
df = @(x) chebSum(APrime_cheb,x);

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
