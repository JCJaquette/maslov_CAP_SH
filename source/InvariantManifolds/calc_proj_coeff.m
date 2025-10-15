function coeff = calc_proj_coeff(eigenvalues, eigenvectors,params)
% This function recursively calculates the coefficients of the
% parameterization for the invariant manifold 
% Inputs: order - order of the desired parameterization 
%         eigenvalues - either the stable or unstable pair of eigenvalues
%         of Df(0)
%         eigenvectors - corresponding eigenvectors
%         params - 2x1 vector in the form [mu, nu']
% Output  coeff - multidimensional array of size (order+1) x (order+1) x 4
%         entry (i,j,k) corresponds to the coefficient k_{i-1,j-1} where in
%         the write up, k is either a,b,c,d and we subtract 1 from each
%         index to be consistent with the zero indexing in the text
%         Each matrix (:,:,i) will be upper triangular (upper left)
    order=params.mfld.order;

    coeff=zeros(order+1,order+1,4);
    if params.isIntval
        coeff = intval(coeff );
    end

    e1=eigenvectors(:,1);
    e2=eigenvectors(:,2);
    lam1=eigenvalues(1);
    lam2=eigenvalues(2);

    
    
    Df0=JacSH(0,params); 
    % the zeroth order coefficient is the equilibrium, corresponding to 0.
    coeff(2,1,:)=e1;
    coeff(1,2,:)=e2;
    
    % suborder=m+n and corresponds to the (m,n)th coefficient
   suborder=2;
   while suborder < order + 1 
        % disp('Calculating terms of order:')
        % disp(suborder)
        for j=suborder:-1:0
            i=suborder-j;
            Aij=(Df0-(i*lam1+j*lam2)*eye(4))^(-1);
            a=coeff(1:i+1,1:j+1,1);
            staraaij=starhat(a,a,i,j);
            staraaaij=tripstarhat(a,a,a,i,j);
            B=[0;0;0;-params.nu*staraaij+staraaaij];
            pij=Aij*B;
            coeff(i+1,j+1,:)=pij;
        end
        suborder=suborder+1;
   end
end