function [count,flag] = getZeroCount(domain_interval, f, f_error, df, df_error, tol, BOOL_plot)

if BOOL_plot

    figure
    hold on
    z = linspace(domain_interval.inf,domain_interval.sup,500);
    for k = 1:500
        fz(k) = real(f(z(k)));
    end
    plot(z,fz,'Color','black')
    hold on

end

count = 0; flag = intval(1)*[];

F = @(x) f(x) + f_error;
dF = @(x) df(x) + df_error;

while isempty(domain_interval) == 0

    b = domain_interval(end);

    if abs(F(b)) > 0

        domain_interval(end) = [];

        if BOOL_plot
            plot([b.inf,b.sup],[0,0],'color','blue')
        end

    else
        
        if F(b.inf)*F(b.sup) < 0 && abs(dF(b)) > 0

            domain_interval(end) = [];
            count = count+1;
            
            if BOOL_plot
                plot([b.inf,b.sup],[0,0],'color','green')
            end

        else

            domain_interval(end) = [];
            end1 = infsup(b.inf,b.mid);
            end2 = infsup(b.mid,b.sup);
            domain_interval = [domain_interval, end1,end2];

        end

    end

    if rad(b)<tol
        domain_interval(end) = [];
        flag = [flag;b];
    end

end



if isempty(flag) == 0
    disp('Not validated, check flag')
    return
end

end

