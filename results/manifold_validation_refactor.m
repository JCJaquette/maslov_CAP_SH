clear all 

%-------------------------------------------------------------------------
% first we calculate the size of the eigenvector/eigenvalue enclosures
% !!! TODO: Make Rigorous; Maybe use formula for e-vec?
% 
point=[0,0,0,0];
MuNu=[0.05,1.6]; % in the order [mu,nu']



params.mu = 0.05; 
params.nu = 1.6; 
params.scale = 1/7;

order = 15-1; 
params.order = order; 
params.mfld.order = order; 

tau = params.scale ;

% Add This boolean for intervals
params.isIntval = 1;


Df0=JacSH(0,MuNu(1),MuNu(2));

[V1,D1]=eigs(Df0); 
[d, ind]=sort(real(diag(D1)));
D=D1(ind,ind);
V=V1(:,ind);

rstar=1e-15;

error=zeros(1,4);
for i=1:4
    error(i)=eig_enclosure(point,MuNu,D(i,i),V(:,i),rstar);
end
error=max(error);

%----------------------------------------------------------------------
% Now we scale the eigenvectors in accordance with Theorem 10.5.1 so the 
% last component is on the order of machine precision.
% TODO: Update Reference
Vscale=V*tau;

% Separate the eigenvalues and associated vectors into those with positive
% and negative real part
uneigs=[D(1,1),D(2,2)];
unvec=Vscale(:,1:2);
stabeigs=[D(3,3),D(4,4)];
stabvec=Vscale(:,3:4);

% -----------------------------------------------------------------------
% Now we calculate the coefficients of the parameterization for the stable
% and unstable manifold up to a desired order. 
order=20;

% unstable
disp('Calculating the coefficients for the unstable manifold.')
uncoeff=calc_proj_coeff(uneigs,unvec,params);
% return

% stable 
disp('Calculating the coefficients for the stable manifold.')
stabcoeff=calc_proj_coeff(stabeigs,stabvec,params);

k=6;
mat=zeros(k^2,4)*stabcoeff(1,1,1);
K = zeros(k^2, 2);
for i = 1:k 
    for j = 1:k
        row=(i-1)*k + j;
        mat(row,:)=stabcoeff(i,j,:);
        K(row,:)=[i-1,j-1];
    end  
end

figure
tiledlayout(2,2) 

nexttile
plot_coeff(uncoeff,order);
title('Coefficient norms for unstable parameterization')

nexttile
plot_coeff(stabcoeff,order);
title('Coefficient norms for stable parameterization')
 
nexttile
plot_manifold(uncoeff,order,'r');
title('Unstable Manifold')

nexttile
plot_manifold(stabcoeff,order,'b');
title('Stable Manifold')


return 
%--------------------------------------------------------------------------
% Now we apply Lemma 4.4 to validate the parameterization we computed 

% TODO: Fix Intvals
 unpoly=radiipoly(MuNu,uncoeff,V,D,order);
 % return
 unstable_bound = min(unpoly(find(unpoly > 0)));

 stpoly=radiipoly(MuNu,stabcoeff,V,D,order);
 stable_bound = min(stpoly(find(stpoly > 0)));
%  
  disp('The error on the parameterization for the unstable manifold is: ')
  disp(unstable_bound)
  disp('The error on the parameterization for the stable manifold is: ')
 disp(stable_bound)

%%%%%%%%%%%%%%%%%%%%%% FUNCTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% This function validates the parameterization of the invariant manifold 
% Inputs: params - 2x1 vector in the form [mu, nu']
%         coeff - multidimensional array of size 3*order+1 x 3*order+1 x 4.
%         This could correspond to the stable or unstable manifold
%         Q - matrix of unscaled eigenvectors 
%         Lambda - diagonal matrix of the eigenvalues 
%         order - order of the approximation - N=m+n
%         rad - the radius on which to search for a negative value
% Output: poly - value at which the sup of the polynomial evaluated with
%         interval arithmetric is negative
function val=radiipoly(params,coeff, Q,Lambda,order)
    % Calculate the coefficient K 
    maxKmn=((order+1)*abs(real(Lambda(1))) - abs(Lambda(1)))^(-1);
 
    K_N = max(abs(Q),[],'all')*max(abs(Q^(-1)),[],'all')*maxKmn;
    
    % extract the coefficients a_{mn}
    
    a=coeff(1:order+1,1:order+1,1);
    
    A=zeros(order+1, 2*order);
    B=zeros(2*order, order+1);
    C=zeros(2*order, 2*order);
    
    bara=[a, A ; B, C];
    %bara=coeff(1:3*order+1, 1:3*order+1,1);
    
    disp('Calculating Y0.')

    %%%%%%
    % Y0 %
    %%%%%%

    quadsum=0;
    suborder=order+1;
    while suborder < 2*order+1
        for i=0:suborder
            j=suborder-i;
            astarij=starhat(bara, bara, i,j);
            quadsum=quadsum+abs(astarij);
        end
        suborder=suborder+1;
    end
   
    cubsum=0;
    suborder=order+1;
    while suborder<3*order+1
        for i=0:suborder
            j=suborder-i;
            tripstara=tripstarhat(bara,bara,bara,i,j);
            cubsum=cubsum+abs(tripstara); 
        end
        suborder=suborder+1;
    end
    
    Y0 = K_N*(params(2)*quadsum+cubsum);
    disp(Y0)
    disp('Calculating Z1.')
    
    %%%%%%
    % Z1 %
    %%%%%%
    
    lowerquadsum=0;
    suborder=1;
    while suborder < order+1
        for i=0:suborder
            j=suborder-i;
            astarij=starhat(bara, bara, i,j);
            lowerquadsum=lowerquadsum+abs(astarij);
        end
        suborder=suborder+1;
    end
    
    linsum=0;
    suborder=1;
    while suborder<order+1
        for i=0:suborder
            j=suborder-i;
            linsum=linsum+abs(bara(i+1,j+1));
        end
        suborder=suborder+1;
    end
    
    Z1= K_N*(2*params(2)*linsum + 3*lowerquadsum);
    disp(Z1);
    disp('Calculating Z2.')
    %%%%%%
    % Z2 %
    %%%%%%
    
    Z2=@(r)K_N*(6*linsum+2*params(2)+3*r);
    disp(Z2)
    
    b=K_N*3;
    a=K_N*(6*linsum + 2*params(2));
  
    disp(a)
    disp(b)
    % this section of code builds an interval on which the sup of the polynomial is negative 
    poly=@(r)(Z1+Z2(r)*r)*r+Y0-r;
    
    %figure(3)
    %fplot(@(r) poly(r));
    
    p=[b, a, Z1-1, Y0];
    
    val=roots(p);
end


