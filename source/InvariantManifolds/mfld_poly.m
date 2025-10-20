function [mflds,r_min] =mfld_poly(params,mflds,BOOL_stable)
% This function validates the parameterization of the invariant manifold 
% Inputs: params - structure
%         mflds  - structure 
%         BOOL_stable - whether to compute the stable manifold (1) or unstable (0). 

% Output: poly - value at which the sup of the polynomial evaluated with
%         interval arithmetric is negative
% 
% 
% OLD
%         Q - matrix of unscaled eigenvectors 
%         Lambda - diagonal matrix of the eigenvalues 
%         order - order of the approximation - N=m+n
%         rad - the radius on which to search for a negative value
%         coeff - multidimensional array of size 3*order+1 x 3*order+1 x 4.
                    %          NOTE-- The dimension of this array seems off
    order=params.mfld.order;

    Q = [mflds.vectors.s,mflds.vectors.u];
    Lambda = diag([mflds.values.s,mflds.values.u]);

    if BOOL_stable
        coeff = mflds.stable.coeffs;
    else
        coeff = mflds.unstable.coeffs;
    end
     

    % Calculate the coefficient K 
    maxKmn= (  (order) *  abs(  real(mflds.values.u(1))  )   );

    K_N = max(abs(Q),[],'all')*max(abs(Q^(-1)),[],'all')/maxKmn;
    
    % extract the coefficients a_{mn}
    
    a=coeff(1:order+1,1:order+1,1);
    
    B=zeros(order+1, 2*order);
    C=zeros(2*order, order+1);
    D=zeros(2*order, 2*order);
    
    bara=[a, B ; C, D];
    %bara=coeff(1:3*order+1, 1:3*order+1,1);

    % Pi_infty 

    Pi_infty=ones(3*order+1, 3*order+1);
    for i = 0:order
        Pi_infty(i+1, 1:(order-i+1) ) = 0*Pi_infty(i+1, 1:(order-i+1) );
    end
    
    disp('Calculating Y0.')

    %%%%%%
    % Y0 %
    %%%%%%

    % Compute a*a

    a_squared =0*bara;

    suborder=1;
    while suborder < 2*order+1
        for i=0:suborder
            j=suborder-i;
            a_squared(i+1,j+1)=starhat(bara, bara, i,j);
        end
        suborder=suborder+1;
    end

    

    % Compute (a*a)*a
    a_cubed =0*bara;
    % suborder=order+1;
    suborder=1;
    while suborder < 3*order+1
        for i=0:suborder
            j=suborder-i;
            a_cubed (i+1,j+1)=starhat(a_squared, bara, i,j);
        end
        suborder=suborder+1;
    end

    quadsum = sum(abs(Pi_infty.*a_squared),'all');
    cubesum = sum(abs(Pi_infty.*a_cubed),'all');



    % % % % % % % 
    % % quadsum=0;
    % % suborder=order+1;
    % % while suborder < 2*order+1
    % %     for i=0:suborder
    % %         j=suborder-i;
    % %         astarij=starhat(bara, bara, i,j);
    % %         quadsum=quadsum+abs(astarij);
    % %     end
    % %     suborder=suborder+1;
    % % end
    % % 
    % % cubsum=0;
    % % suborder=order+1;
    % % while suborder<3*order+1
    % %     for i=0:suborder
    % %         j=suborder-i;
    % %         tripstara=tripstarhat(bara,bara,bara,i,j);
    % %         cubsum=cubsum+abs(tripstara); 
    % %     end
    % %     suborder=suborder+1;
    % % end
    
    Y0 = K_N*(params.nu*quadsum+cubesum);
    disp(Y0)
    disp('Calculating Z1.')
    
    %%%%%%
    % Z1 %
    %%%%%%
    
    % lowerquadsum=0;
    % suborder=1;
    % while suborder < order+1
    %     for i=0:suborder
    %         j=suborder-i;
    %         astarij=starhat(bara, bara, i,j);
    %         lowerquadsum=lowerquadsum+abs(astarij);
    %     end
    %     suborder=suborder+1;
    % end
    % 
    % linsum=0;
    % suborder=1;
    % while suborder<order+1
    %     for i=0:suborder
    %         j=suborder-i;
    %         linsum=linsum+abs(bara(i+1,j+1));
    %     end
    %     suborder=suborder+1;
    % end
    
    linsum       = sum( abs(bara .* ( ~Pi_infty )) , 'all');
    lowerquadsum = sum( abs(a_squared .* ( ~Pi_infty )) , 'all');

    
    Z1= K_N*(2*params.nu*linsum + 3*lowerquadsum);
    disp(Z1);
    disp('Calculating Z2.')
    %%%%%%
    % Z2 %
    %%%%%%
    
    Z2=@(r)K_N*(6*linsum+2*params.nu+3*r);
    disp(Z2)
    
    b=K_N*3;
    a=K_N*(6*linsum + 2*params.nu);
  
    disp(a)
    disp(b)
    % this section of code builds an interval on which the sup of the polynomial is negative 
    poly=@(r)(Z1+Z2(r)*r)*r+Y0-r;
    
    %figure(3)
    %fplot(@(r) poly(r));
    
    % get numerical roots
    if params.isIntval
        p=[b.sup, a.sup, Z1.sup-1, Y0.sup];
    else
        p=[b, a, Z1-1, Y0];
    end
    
    cube_roots =roots(p);


    r_min = min(cube_roots (find(cube_roots > 0)));
    if r_min == 0
        r_min = 10*Y0;
    end

    % Check that there is a positive root
    if isempty(r_min)
        r_min = NaN;
    else
        % inflate a bit
        r_min = r_min*1.00001;
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


end
