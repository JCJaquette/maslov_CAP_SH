function [bndl_local] = bndl_one_point(phi1, phi2, bndl, params)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
coeff_size  = size(bndl.coeffs);


bndl_local=[]; 
for i=1:4
    column_coeff =  bndl.coeffs(:,i,:,:);
    column_coeff = permute(column_coeff ,[3,4,1,2]);

    column_coeff = reshape(column_coeff ,coeff_size(3),coeff_size(4),4);

    val = mfld_one_point(phi1, phi2, column_coeff, params);
    bndl_local= [bndl_local,val];
end


end


