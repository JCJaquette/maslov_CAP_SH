function [mflds,r_min, data_mfld_poly ] =mfld_poly(params,mflds,BOOL_stable)
% This function validates the parameterization of the invariant manifold 
% Inputs: params - structure
%         mflds  - structure 
%         BOOL_stable - whether to compute the stable manifold (1) or unstable (0). 

% Output: mflds - mflds, updated with the validation radius
%         r_min - value for which p(r_min) < 0 ; or NaN 
% 

    order=params.mfld.order;

%   Q - matrix of unscaled eigenvectors 
    Q = [mflds.vectors.s,mflds.vectors.u];

    if BOOL_stable
        coeff = mflds.stable.coeffs;
    else
        coeff = mflds.unstable.coeffs;
    end

    % Calculate the coefficient K_N 
    maxKmn= (  (order) *  abs(  real(mflds.values.u(1))  )   );
    Q_inv = Q^(-1);
    Q_inv_last = abs(Q_inv(:,4));

    K_N = norm(Q,1)*max(Q_inv_last)/maxKmn
    data_mfld_poly.K_N =K_N ;
    
    % extract the coefficients a_{mn}
    a=coeff(1:order+1,1:order+1,1);    
    B=zeros(order+1, 2*order);
    C=zeros(2*order, order+1);
    D=zeros(2*order, 2*order);
    
    % Makes a vector of size 3*N x 3*N.
    bar_a=[a, B ; C, D];
    

    % Pi_infty 
    Pi_infty=ones(3*order+1, 3*order+1);
    for i = 0:order
        Pi_infty(i+1, 1:(order-i+1) ) = 0*Pi_infty(i+1, 1:(order-i+1) );
    end
    
    

    %%%%%%
    % Y0 %
    %%%%%%
    disp('Calculating Y0.')

    % Compute a*a
    a_squared =0*bar_a;

    suborder=1;
    while suborder < 2*order+1
        for i=0:suborder
            j=suborder-i;
            a_squared(i+1,j+1)=starhat(bar_a, bar_a, i,j);
        end
        suborder=suborder+1;
    end

    

    % Compute (a*a)*a
    a_cubed =0*bar_a; 
    suborder=1;
    while suborder < 3*order+1
        for i=0:suborder
            j=suborder-i;
            a_cubed (i+1,j+1)=starhat(a_squared, bar_a, i,j);
        end
        suborder=suborder+1;
    end

    quadsum = sum(abs(Pi_infty.*a_squared),'all');
    cubesum = sum(abs(Pi_infty.*a_cubed),'all');
    
    Y0 = K_N*(params.nu*quadsum+cubesum);
    disp(Y0)
    data_mfld_poly.Y0 =Y0 ;
    
    %%%%%%
    % Z1 %
    %%%%%%
    disp('Calculating Z1.')

    linsum       = sum( abs(bar_a ) , 'all');
    lowerquadsum = sum( abs(a_squared ) , 'all');

    Z1= K_N*(2*params.nu*linsum + 3*lowerquadsum);
    disp(Z1);
    data_mfld_poly.Z1=Z1;

    
    %%%%%%
    % Z2 %
    %%%%%%
    disp('Calculating Z2.')
    Z2_0=K_N*(6*linsum + 2*params.nu);
    Z2_1=K_N*3;
    
    Z2=@(r)K_N*(6*linsum+2*params.nu+3*r);
    disp(Z2)
    
  
    disp(Z2_0)
    disp(Z2_1)

    data_mfld_poly.Z2_0 = Z2_0;
    data_mfld_poly.Z2_1 = Z2_1;

    % this section of code builds an interval on which the sup of the polynomial is negative 
    poly=@(r)(Z1+Z2(r)*r)*r+Y0-r;
    
    data_mfld_poly.poly = poly;
    
    % get numerical roots
    if params.isIntval
        p=[Z2_1.sup, Z2_0.sup, Z1.sup-1, Y0.sup];
    else
        p=[Z2_1, Z2_0, Z1-1, Y0];
    end
    
    cube_roots =roots(p);


    r_min = min(cube_roots (find(cube_roots >= 0)));
    if r_min == 0
        r_min = 10*Y0;
    end
    

    % Check that there is a positive root
    if isempty(r_min)
        r_min = NaN;
    else
        % inflate a bit
        r_min = r_min*1.00001;

        if params.isIntval
            r_min  = intval(r_min);
        end

        % ensure the root is indeed negative
        if ~(poly(r_min)<0)
            r_min=NaN;
        end
    end

    % Define the validation radius (Default is NaN)
    if BOOL_stable
        mflds.stable.r_min = r_min;
    else
        mflds.unstable.r_min = r_min;
    end

    data_mfld_poly.r_min = r_min;

end
