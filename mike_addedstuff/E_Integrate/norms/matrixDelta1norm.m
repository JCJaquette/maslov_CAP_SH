function normout = matrixDelta1norm(mat,del)
% matrix norm for weighted 1-norm

    [~,n] = size(mat);

    normout = 0;

    for i = 1:n        
        if vectorDelta1norm(mat(:,i),del) /del^i > normout
            normout = vectorDelta1norm(mat(:,i),del) /del^i;
        end
    end

end

