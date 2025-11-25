function [D1, D2] = differentiate_mani(coeffs)

    D1 = coeffs;
    D2 = coeffs;
    ord = length(coeffs(:,1,1));

    for i = 1:ord

        D1(i,:,:) = (i-1)*D1(i,:,:);
        D2(:,i,:) = (i-1)*D2(:,i,:);

    end

    D1 = [D1(2:end,:,:);0*D1(1,:,:)];
    D2 = [D2(:,2:end,:),0*D2(:,1,:)];

end

