function [out] = avgnorms(mani)

    normsu = zeros(1,59);
    for i = 1:59
    c1 = reshape(mani(60,i,:),[4,1]);
    c2 = reshape(mani(60,i+1,:),[4,1]);
    normsu(i) = norm(c1-c2);
    end

    out = mean(normsu);

end

