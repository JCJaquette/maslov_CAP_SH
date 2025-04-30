function normout = matrixDelta1norm(mat,del)
% matrix norm for weighted 1-norm

    [n,~] = size(mat);

    delvec = del.^(0:1:n-1);

    weightedcolsums = delvec*mat; % Shouldn't we need to take an absolute value of mat here? 

    [mx,ind] = max(weightedcolsums);

    normout = mx/del^(ind-1);

end

