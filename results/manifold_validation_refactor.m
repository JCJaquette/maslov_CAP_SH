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

order = 15-1;  % Manifold
params.order = order; 
params.mfld.order = order; 

tau = params.scale ;

% Add This boolean for intervals
params.isIntval = 0;


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


  
%--------------------------------------------------------------------------
% Now we apply Lemma 4.4 to validate the parameterization we computed 

% TODO: Fix Intvals
 unpoly=mfld_poly(params,uncoeff,V,D);
 % return
 unstable_bound = min(unpoly(find(unpoly > 0)));

 stpoly=mfld_poly(params,stabcoeff,V,D);
 stable_bound = min(stpoly(find(stpoly > 0)));
%  
  disp('The error on the parameterization for the unstable manifold is: ')
  disp(unstable_bound)
  disp('The error on the parameterization for the stable manifold is: ')
 disp(stable_bound)
