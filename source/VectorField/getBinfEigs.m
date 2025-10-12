function [vectors, values] = getBinfEigs(params)

    % % % mat=B_infinity(params);
    % % % try 
    % % %     [V1,D1]=eigs(mat);
    % % % catch
    % % %     [V1,D1]=eigs(mid(mat));
    % % % end
    % % % % sort the eigenvalues and eigenvectors by increasing real part
    % % % [d, ind]=sort(real(diag(D1)));
    % % % D=D1(ind,ind);
    % % % V=V1(:,ind);
    % % % 
    % % % % scale the eigenvectors if necessary
    % % % 
    % % % 
    % % % % save the eigenvalues and eigenvectors
    % % % values.s=[D(1,1),D(2,2)];
    % % % vectors.s=V(:,1:2);
    % % % values.u=[D(3,3),D(4,4)];
    % % % vectors.u=V(:,3:4);
    
    % % % if ((real(values.s(1))<0) == 1) == 0
    % % %     disp('There may be something wrong with the assignment of the asymptotic stable and unstable vectors')
    % % % end

    % New/explicit  Approach 
    [rho theta]= rTh_coord(params) ;
    r=rho;
    
    mu_s1 = - sqrt(r)*exp(1i*theta/2);
    mu_u1 =   sqrt(r)*exp(1i*theta/2);
    mu_s2 = - sqrt(r)*exp(-1i*theta/2);
    mu_u2 =   sqrt(r)*exp(-1i*theta/2);

% Define eigenvector: cf (3.6) 
    Vu1 = [exp(-1i*theta) / r ; 
            1 ; 
            sqrt(r) *exp(1i*theta/2)+2*exp(-1i*theta/2)/sqrt(r) ; 
            exp(-1i*theta/2)/sqrt(r) ];

% scale the eigenvectors if necessary
    Vu1 = Vu1/mid(sqrt( 1+1/r^2+5/r+r+4*cos(theta) ) );

    Vu2 = conj(Vu1);
    Vs1 = Vu1.*[1;1;-1;-1];
    Vs2 = conj(Vs1);

% Store data
    values.s = [mu_s1,mu_s2];
    vectors.s = [Vs1,Vs2];

    values.u = [mu_u1,mu_u2];
    vectors.u = [Vu1,Vu2];


        
end