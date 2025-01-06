function normout = vectorDelta1norm(vec,del)
% weighted 1-norm

    N = length(vec);

    for i = 1:N
        vec(i) = del^(i-1)*vec(i);
    end

    normout = norm(vec,1);

end

