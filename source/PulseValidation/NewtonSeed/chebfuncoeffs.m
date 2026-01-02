function [x] = chebfuncoeffs(pSoln,order)
    
    a1 = chebfun(pSoln(:,1),'equi');
    x.a1 = a1.funs{1,1}.onefun.coeffs;
    a2 = chebfun(pSoln(:,2),'equi');
    x.a2 = a2.funs{1,1}.onefun.coeffs;
    a3 = chebfun(pSoln(:,3),'equi');
    x.a3 = a3.funs{1,1}.onefun.coeffs;
    a4 = chebfun(pSoln(:,4),'equi');
    x.a4 = a4.funs{1,1}.onefun.coeffs;   

    x.a1 = [x.a1(1:min(length(x.a1),order));zeros(order - length(x.a1),1)]';
    x.a2 = [x.a2(1:min(length(x.a2),order));zeros(order - length(x.a2),1)]';
    x.a3 = [x.a3(1:min(length(x.a3),order));zeros(order - length(x.a3),1)]';
    x.a4 = [x.a4(1:min(length(x.a4),order));zeros(order - length(x.a4),1)]';
    x.a1(2:end) = x.a1(2:end)/2;
    x.a2(2:end) = x.a2(2:end)/2;
    x.a3(2:end) = x.a3(2:end)/2;
    x.a4(2:end) = x.a4(2:end)/2;

end

