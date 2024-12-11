function normout = matrixDelta1norm(mat,del)

    [~,n] = size(mat);

    normout = 0;

    for i = 1:n        
        if vectorDelta1norm(mat(:,i),del) /del^i > normout
            normout = vectorDelta1norm(mat(:,i),del) /del^i;
        end
    end

end

