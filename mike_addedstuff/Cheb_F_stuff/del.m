function [out] = del(i,j)
% delta function

    if i ~= j
        out = 0;
    else
        out = 1;
    end

end

