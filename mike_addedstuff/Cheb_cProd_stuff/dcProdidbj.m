function [out] = dcProdidbj(a,i,j)
% doing DcPProd the same way as was done with DcProd was too slow
% dcp is the same fn as chebcPProd_kth, but is a bit faster by not summing things times 0

    b = 0*a;

    b(j) = 1;

    out = dcp(a',a',b',i);

end
