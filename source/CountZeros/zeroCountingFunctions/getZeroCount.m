function [count,flag] = getZeroCount(domain_interval, f, f_error, df, df_error, tol, BOOL_plot)

if BOOL_plot
%plotting the function
    if isintval(f(0))
        figure
        hold on
        z = linspace(domain_interval.inf,domain_interval.sup,500);
        for k = 1:500
            fzk = f(z(k));
            fz(k) = mid(real(fzk));
        end
        plot(z,fz,'Color','black')
        hold on
    else
        figure
        hold on
        z = linspace(domain_interval.inf,domain_interval.sup,500);
        for k = 1:500
            fz(k) = real(f(z(k)));
        end
        plot(z,fz,'Color','black')
        hold on
    end

end

%counting zeros 
count = 0; flag = intval(1)*[];

F = @(x) f(x) + f_error;
dF = @(x) df(x) + df_error;

list_domain = intval(1)*[];
list_f = intval(1)*[];

while isempty(domain_interval) == 0

    b = domain_interval(end);

    Fb = F(b);
    if abs(Fb) > 0
        
        list_domain(end+1) = domain_interval(end);
        list_f(end+1) = Fb;
        domain_interval(end) = [];

        if BOOL_plot
            plot([b.inf,b.sup],[0,0],'color','blue')
        end

    else
        
        dfb = dF(b); 
        endpt_prod = F(b.inf)*F(b.sup);

        if endpt_prod < 0 && abs(dfb) > 0

            list_domain(end+1) = domain_interval(end);
            list_f(end+1) = Fb;
            domain_interval(end) = [];
            count = count+1;
            
            if BOOL_plot
                plot([b.inf,b.sup],[0,0],'color','green')
            end

        elseif endpt_prod > 0 && abs(dfb) > 0

            list_domain(end+1) = domain_interval(end);
            list_f(end+1) = Fb;
            domain_interval(end) = [];
    
            if BOOL_plot
                plot([b.inf,b.sup],[0,0],'color','blue')
            end

        else

            if rad(b)<tol
                domain_interval(end) = [];
                flag = [flag;b];
            else
                domain_interval(end) = [];
                end1 = infsup(b.inf,b.mid);
                end2 = infsup(b.mid,b.sup);
                domain_interval = [domain_interval, end1,end2];

            end

        end

    end

end

if BOOL_plot
    figure
    plot(list_domain,real(list_f))

end

if isempty(flag) == 0
    disp('Not validated, check flag')
    
    return
end

end