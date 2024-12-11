function normout = vectorDelta1norm(vec,del)

    N = length(vec);

    for i = 1:N
        vec(i) = del^(i-1)*vec(i);
    end

    normout = norm(vec,1);

end

