function [mpoint,min,mps] = closestpt(mani,point)

    n = length(mani(:,1,1));
    m = length(mani(1,:,1));

    min = 100;
    c = [0;0;0;0];

    for i = 1:n
        for j = 1:m

            for k = 1:4
                c(k) = mani(i,j,k);
            end

            if norm(c-point,2) < min
                min = norm(c-point);
                mpoint = mani(i,j,:);
                mps = [i,j];
            end

        end
    end

    

end

