function I = RadiiPoly1D(f,df,xbar,dom)
% T(x) = x - f(x)/A 
% Bound |T(xbar)|, |DT(xbar)|
    
    A = intval(1)*df(xbar);
    DT = @(x) 1 - df(x)/A;

    Y0 = abs(f(xbar)/A);

    max_potential_r = min(xbar.sup - dom(1), dom(end) - xbar.inf);
    grow = max_potential_r/30;

    DT_bound = 0;
    grow1 = xbar; grow_1 = xbar;
    good_rs = []; r = 0;
    stop = 0;

    while isempty(good_rs) == 1

        for i = 1:30
    
            r = r + grow;
            grow1 = infsup(grow1.sup, grow1.sup + grow);
            grow_1 = infsup(grow_1.inf - grow, grow_1.inf);
    
            DT_on_g1 = abs(DT(grow1));
            DT_on_g_1 = abs(DT(grow_1));
    
            DT_bound = max([DT_bound, DT_on_g1.sup, DT_on_g_1.sup]);
    
            if (DT_bound - 1)*r + Y0 < 0
    
                good_rs = [good_rs,r];
    
            end
    
            if isempty(good_rs) == 0 && (DT_bound - 1)*r + Y0 > 0
                break
            end
    
        end

        grow = grow/2; 
        grow1 = xbar; grow_1 = xbar;

        stop = stop+1;
        if stop == 20
            disp('radii polynomial failed for x = ')
            xbar
            break
        end

    end

    if isempty(good_rs) == 0
        I = infsup(good_rs(1), good_rs(end));
    else
        I = intval(0);
    end
    
end

