function [bndl]= getAllBundleCoefficients(params,mflds)
% TODO: Restructure this so that the manifold coefficients are already
% computed. 
% TODO: Allow for either stable or unstable manifold to be called
    order = params.mfld.order; 
 %   eq = zeros(1, 4); 
 %   Df0 = JacSH_toMerge(0, params); 


 % Add Coeff to manifold data structure
values = mflds.values;
vectors = mflds.vectors;


    


    %%% NOTE! See ordering of eigenvalues, (3.5)
    e_val = [values.s, values.u];

    disp('Computing Manifold')
    G_hat = DFQbundle(params, mflds); % I think this is \hat G

    %%% NOTE! See ordering of eigenvalues, (3.5)
    Q0 = [vectors.s, vectors.u]; 

    % Rescale the eigenvectors as needed
    Q0 = Q0*diag(1./Q0(2,:)); %%%% Why this??
    
    % original coordinates 
    coeffs = zeros(4,4,order + 1, order + 1); % W--series
    if params.isIntval
        coeffs =intval(coeffs);
    end
    coeffs(:,:,1,1)=Q0;

    

    coeffs_tilde = coeffs; % W tilde --series
    coeffs_tilde(:,:,1,1)=eye;

    
    normalForm_coeff = zeros(4,4,order + 1, order + 1); % A
    % normalForm_coeff(:,:,1,1)=Omega;

    if params.isIntval
        coeffs =intval(coeffs);
        normalForm_coeff=intval(normalForm_coeff);
    end
    
    disp('Computing Bundles')
    for alpha = 1:order
        % alpha
        for n = 0:alpha 
            
            
            m = alpha - n;


            %s_ij = starhatMat(A,coeffs, i,j); 
            % disp(['(',num2str(m),',',num2str(n),')'])
            % disp(j)
            
            % equation 6.13
            % Since the (m,n) component of W is zero, this is
            % equivalent to 
            %       G_\alpha . W_0 + G \hat * W
            s_mn = starMat(G_hat,coeffs, m,n); 
            % Add the normal form part
            s_mn=s_mn-starhatMat(coeffs,normalForm_coeff, m,n);  

            % 4 x 4 vector
            s_tilde_mn = Q0\s_mn; % S tilde

            % known resonances
            if alpha ==2 
                % keyboard
                if n==0 %(2,0)
                    normalForm_coeff(1,3,m+1,n+1) = s_tilde_mn(1,3);
                elseif n==1 %(1,1)
                    normalForm_coeff(1,4,m+1,n+1) = s_tilde_mn(1,4);
                    normalForm_coeff(2,3,m+1,n+1) = s_tilde_mn(2,3);
                elseif n==2 %(0,2)
                    normalForm_coeff(2,4,m+1,n+1) = s_tilde_mn(2,4);
                end
            end

            coeff_new = K_op(e_val,m,n).*s_tilde_mn; 

            % eigenbasis_coeffs(:, :, i + 1, j + 1) = [coeff1, coeff2];
            coeffs(:,:,m+1, n+1) = Q0*coeff_new; 

        end
    end

    bndl.coeffs = coeffs;
    bndl.normalForm = normalForm_coeff;


end



