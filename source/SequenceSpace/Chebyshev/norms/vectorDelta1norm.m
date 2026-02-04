function normout = vectorDelta1norm(vec,del)
% weighted 1-norm

    N = length(vec);
    [m,~] = size(vec);

    if m == 1
        weightedvec = (del.^(0:N-1)).*vec;
    else
        weightedvec = (del.^(0:N-1)').*vec;   
    end

    normout = norm(weightedvec,1);

end

