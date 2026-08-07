function flag = intersectioncheck(dom,F,dF,a1,a4,b1,b4)

flag = 'Problem occurred with hypothesis 3.17';

while ~isempty(dom)

    b = dom(end);
    endpt_prod = F(b.inf)*F(b.sup);

    if endpt_prod < 0 && abs(dF(b)) > 0

        if abs(a1(b)) > 0 || abs(a4(b)) > 0 || abs(b1(b)) > 0 || abs(b4(b)) > 0

            dom(end) = [];

        else

                dom(end) = [];
                end1 = infsup(b.inf,b.mid);
                end2 = infsup(b.mid,b.sup);
                dom = [dom, end1,end2];

        end

    else

        if abs(F(b)) > 0

            dom(end) = [];

        end
        
    end

end

if isempty(dom)
    flag = 'Hypothesis 3.17 verified';
end

end