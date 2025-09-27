function [coeffs, normalForm_coeff,mflds]= getAllBundleCoefficients(params)
% TODO: Restructure this so that the manifold coefficients are already
% computed. 
    order = params.mfld.order; 
 %   eq = zeros(1, 4); 
 %   Df0 = JacSH_toMerge(0, params); 

 % Add Coeff to manifold data structure

    vectors.s = params.eigenvectors.s; 
    values.s = params.eigenvalues.s;
    vectors.u = params.eigenvectors.u; 
    values.u = params.eigenvalues.u;

   % v0 = vectors.s(:, 1)*params.scale; 
    lam1 = values.s(1); 
    lam2 = values.s(2); 

    % lam2u = values.u(1); 
    % lam1u = values.u(2); 

    % e_val = [lam1 lam2 lam1u lam2u];

    e_val = [lam1 lam2 fliplr(values.u)];%%% NOTE! See ordering of eigenvalues, (3.5)

    
    mflds.coeff.s = calc_proj_coeff(values.s, vectors.s, params); 
    G_hat = DFQbundle(params, mflds); % I think this is \hat G


    Q0 = [vectors.s, fliplr(vectors.u)]; %%% NOTE! See ordering of eigenvalues, (3.5)

    % Rescale the eigenvectors as needed
    Q0 = Q0*diag(1./Q0(2,:)); %%%% Why this??
    
    % original coordinates 
    coeffs = zeros(4,4,order + 1, order + 1); % W--series
    coeffs(:,:,1,1)=Q0;

    coeffs_tilde = coeffs; % W tilde --series
    coeffs_tilde(:,:,1,1)=eye;

    
    normalForm_coeff = zeros(4,4,order + 1, order + 1); % A
    % normalForm_coeff(:,:,1,1)=Omega;

    
    for alpha = 1:order
        for n = 0:alpha 
            
            
            m = alpha - n;


            %s_ij = starhatMat(A,coeffs, i,j); 
            disp(['(',num2str(m),',',num2str(n),')'])
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

    Za =  sum(abs(G_hat),'all')/order
    Zc =  sum(abs(normalForm_coeff),'all')/order
end


function K_alpha = K_op(e_val,m,n)
% K_alpha - but only multiplied against the coefficients, not the normal
% form. 
% Resonances are set to zero. 
    % e_val -- a row vector
    alpha = m+n;

    denominator = m*e_val(1)+n*e_val(2)  +   e_val-transpose(e_val);
    K_alpha = 1./denominator ;

    if alpha ==2 
        if n==0 %(2,0)
            K_alpha(1,3)=0;
        elseif n==1 %(1,1)
            K_alpha(1,4)=0;
            K_alpha(2,3)=0;
        elseif n==2 %(0,2)
            K_alpha(2,4)=0;
        end
    end
end
