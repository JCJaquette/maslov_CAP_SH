% x is a 4 dim vector
function [vectors, values]= getJacEigs_explicit(params)

% This is how we compute the eigenvectors, although we use a simplification
% to minimize the error a bit
        % % % Define change of coordinates  symplectic = S * normal
        % % S= [1 0 0 0;
        % %     0 0 1 0;
        % %     0 2 0 1;
        % %     0 1 0 0];
        % % S_inv...
        % %  = [1 0 0 0;
        % %     0 0 0 1;
        % %     0 1 0 0;
        % %     0 0 1 -2];
        % % % Get explicit eigenvalues from B_inf defined in (3.2)
        % % [vectors_new, values_new] = getBinfEigs(params);
        % % vectors_new.u = S_inv*vectors_new.u;
        % % vectors_new.s = S_inv*vectors_new.s;

   
% Get r theta coordinates
    [rho, theta]= rTh_coord(params) ;
    r=rho;
% Define eigenvalues
    mu_s1 = - sqrt(r)*exp(1i*theta/2);
    mu_u1 =   sqrt(r)*exp(1i*theta/2);
    mu_s2 = - sqrt(r)*exp(-1i*theta/2);
    mu_u2 =   sqrt(r)*exp(-1i*theta/2);

% Define eigenvector: cf (3.6), although multiplied by S^-1
    Vu1 = [exp(-1i*theta) / r ; 
            exp(-1i*theta/2)/sqrt(r) ;
            1 ; 
            sqrt(r) *exp(1i*theta/2) ];

 % scale the eigenvectors  
     if params.isIntval
        norm_vec = mid(norm(Vu1));
     else 
         norm_vec = norm(Vu1);
     end
    Vu1 = Vu1/norm_vec;
    % Vu1 =*Vu1;

    Vu2 = conj(Vu1);
    Vs1 = Vu1.*[1;-1;1;-1];
    Vs2 = conj(Vs1);
% Store data
    vectors.u = [Vu1 , Vu2 ] ;
    vectors.s = [Vs1 , Vs2 ] ;

    values.u = [mu_u1, mu_u2];
    values.s = [mu_s1, mu_s2];
    
end
