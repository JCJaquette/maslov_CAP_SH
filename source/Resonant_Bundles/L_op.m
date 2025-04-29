function [ L_entry,L_mat ] = L_op(N,mu_s,mu)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here
% Computes K^s a_in
% mu_s - eigenvalue

% Leading Operator
z_type = 0*mu_s +0*mu; % Zero, but potentially an interval

L_entry = z_type *zeros(N+1,N+1);
for m=0:N
    for n=0:N
        xi = -(1+((m-1)*mu_s+n*mu_s')^2)^2-mu;
        if (m==2&&n==0)||(m==1&&n==1)
            xi=z_type; % Going to be zero anywhay
        elseif m==1&&n==0
            xi = -1-(mu+z_type);
        end
        if isnan(xi)
            keyboard
        end
        L_entry(m+1,n+1)=xi;
    end
end

L_mat=reshape(L_entry,(N+1)^2,1); 


L_mat=diag(L_mat); 

end