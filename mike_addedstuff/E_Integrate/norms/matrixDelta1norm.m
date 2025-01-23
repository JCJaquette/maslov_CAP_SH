function normout = matrixDelta1norm(mat,del)
% matrix norm for weighted 1-norm

    [n,~] = size(mat);

    delvec = del.^(0:1:n-1);

    weightedcolsums = delvec*mat;

    [mx,ind] = max(weightedcolsums);

    normout = mx/del^(ind-1);

end

