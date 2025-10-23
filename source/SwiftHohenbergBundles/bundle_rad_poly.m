function [ r_min, data_bndl_poly] = bundle_rad_poly(params,mflds,bndl)
% This function validates the parameterization of the bundles of the invariant manifold 
% Inputs: params - structure
%         mflds  - structure 
%         BOOL_stable - whether to compute the stable manifold (1) or unstable (0). 
% Output: 
%         r_min - value for which p(r_min) < 0 ; or NaN 
%         bndl_poly - contains all the data from the rad poly calculation 
% 

    % Intlab compatible
    if params.isIntval
        zero=intval(0);
    else
        zero=0;
    end


%   Instantiate several variables 
    All_Bundle_coeffs = bndl.coeffs;
    normalForm_coeff = bndl.normalForm;
    order = params.mfld.order; 
    values = mflds.values;
    e_val = [values.s, values.u];

    % When products are taken, everything will be inside this order.
    order_2 = 3*(order+1)-1;

    params_extend = params;
    params_extend.mfld.order=order_2 ;
    mflds_extend =mflds;

    % Extend things for taking products 
    P =  zeros(order_2 +1,order_2 +1,4)*zero;
    P(1:order +1,1:order +1,:)=mflds.stable.coeffs;
    mflds_extend.stable.coeffs=P;

    All_Bundle_coeffs_extend =  zeros(4,4,order_2 +1,order_2 +1)*zero;
    All_Bundle_coeffs_extend(:,:,1:order +1,1:order +1) =All_Bundle_coeffs;
 
    % W_0 -- The zeroth bundle coefficient
    W_0 = All_Bundle_coeffs(:,:,1,1);
    W_0_inv = inv(W_0);

    W_0_norm = norm(W_0 ,1);
    W_0_inv_norm = norm(W_0_inv ,1);

    % Calcuate \hat G^N
    G_hat_N = DGPbundle(params_extend, mflds_extend); 
 
 
    
    p1 = mflds.stable.coeffs(:,:,1);


    

    % P_infty
    P_infty = + mflds.stable.r_min;


    % Epsilon_infty
    % Lemma 6.7

    eps_infty = 2*params.nu* P_infty +3*(2*sum(abs(p1),'all')*P_infty +P_infty^2 );
    
    % K_N 
    % Eq (3.7)
    [rho, theta]= rTh_coord(params) ;


    K_N_s = 1 / ( (order+1) * sqrt(rho) * cos( theta/2) );
    K_N_u = 1 / ( order * sqrt(rho) * cos( theta/2) );
    K_Np2_u = 1 / ( (order +2)* sqrt(rho) * cos( theta/2) );

%%% More bounds 
G_N_norm = sum(abs(G_hat_N),'all'); 
W_N_norm = sum(abs(All_Bundle_coeffs),'all'); 

    %%%%%%%%%%%%%%%%%
    %%%  Y bound  %%% 
    %%%%%%%%%%%%%%%%%

    %%%%%%%%%%%%%%%
    %%%  Y_0^a  %%%
    %%%%%%%%%%%%%%%
    disp('Computing Y bound')
    
    little_sum_2 = zero; 

    vault =  zeros(4,3*order +1,3*order +1)*zero;
    for alpha = order+1:3*order 
        for i = 0:alpha 
            j = alpha - i;  
            summand_2 = K_op(e_val,i,j).* (W_0\ starMat(G_hat_N,All_Bundle_coeffs_extend,i,j));
            vault(:,i+1,j+1)=sum(abs(summand_2)) ;
            little_sum_2 =little_sum_2 +norm(abs(summand_2),1); 
        end
    end

    Ya_0_vec =zero*(1:4);
    for j = 1:4
        Ya_0_vec(j)=sum(vault(j,:,:),'all');
    end
    Ya_0 = max(Ya_0_vec);
    try
        Ya_0_vec.sup
        Ya_0_vec.inf
        Ya_0_vec.rad
    end

    data_bndl_poly.Ya_0_vec=Ya_0_vec;

 


    %%%%%%%%%%%%%%%
    %%%  Y_0^b  %%%
    %%%%%%%%%%%%%%%

    Yb_0_vec = (1+zero*Ya_0_vec)  ; % vector of ones 
    Yb_0_vec = Yb_0_vec  * W_0_inv_norm * eps_infty* W_N_norm;

    Yb_0_vec(1:2)=K_N_s*Yb_0_vec(1:2);
    Yb_0_vec(3:4)=K_N_u*Yb_0_vec(3:4);
    

    Yb_0 = max(Yb_0_vec);
    data_bndl_poly.Yb_0_vec=Yb_0_vec;

    %%%%%%%%%%%%%%%
    %%%  Y_0^c  %%%
    %%%%%%%%%%%%%%%
% This only needs things of order N+1 and N+2

    Yc_0_vec = zero*Ya_0_vec  ; % vector of zero
    
normalForm_coeff;

    little_sum = zero*zeros(4,4); 

    % Beta is order 2 
    for i_beta = 0:2 
        j_beta = 2-i_beta ;
        A_beta = abs(normalForm_coeff(:,:,i_beta+1,j_beta+1));

        for alpha = order+1:order+2
            for i = 0:alpha 
                j = alpha - i;  

                if (i+1-i_beta > 0) && (j+1-j_beta > 0)
                    W_a_minus_b = abs(All_Bundle_coeffs(:,:,i+1-i_beta,j+1-j_beta));
                else 
                    continue
                end

                summand = abs(K_op(e_val,i,j)).* (W_a_minus_b * A_beta);
                little_sum =little_sum + summand ;
            end
        end

    end
    Yc_0_vec = sum(little_sum );
    data_bndl_poly.Yc_0_vec=Yc_0_vec;
    
    Yc_0 = max(Yc_0_vec )

 
    %%%%%%%%%%%%%%%%%
    %%%  Z bound  %%% 
    %%%%%%%%%%%%%%%%%

    %%%%%%%%%%%%%%%
    %%%  Z_0^a  %%% 
    %%%%%%%%%%%%%%%

    % first_component = sum(abs(G_hat_N(:,:,1,1)),'all');

    vault =  zeros(4,order_2+1,order_2+1)*zero;

    little_sum = zero; 
    for alpha = 1:order_2 
        % K_N_alpha = 1 / (( order+alpha) * sqrt(rho) * cos( theta/2) ) ;
        for i = 0:alpha 
            j = alpha - i;  
            % summand = starMat(G_hat_N,All_Bundle_coeffs_extend,i,j);
            summand =  (W_0\ G_hat_N(:,:,i+1,j+1)*W_0 );
            vault(:,i+1,j+1)= sum(abs(summand));
            % little_sum =little_sum +sum(abs(summand),'all'); 
            little_sum =little_sum +norm(abs(summand),1); 
            % little_sum_2 =little_sum_2 +norm(abs(summand),1)*K_N_alpha ; 
        end
    end

    Za_vec =zero*(1:4);
    for j = 1:4
        Za_vec (j)=sum(vault(j,:,:),'all');
        if j <= 2 
            Za_vec (j) = Za_vec(j) * K_N_s;
        else
            Za_vec (j) = Za_vec(j) * K_N_u;
        end
    end
    data_bndl_poly.Za_vec =Za_vec ;
    Za = max(Za_vec )

    %%%%%%%%%%%%%%%
    %%%  Z_0^b  %%%
    %%%%%%%%%%%%%%%

    Zb_vec = (1+zero)*Za_vec  ; % vector of ones 
    Zb_vec  = Zb_vec   * eps_infty*W_0_inv_norm*W_0_norm;

    Zb_vec(1:2)=K_N_s*Zb_vec(1:2);
    Zb_vec(3:4)=K_N_u*Zb_vec(3:4);

    Zb = max(Zb_vec );
    data_bndl_poly.Zb_vec =Zb_vec ;

    %%%%%%%%%%%%%%%
    %%%  Z_0^c  %%% 
    %%%%%%%%%%%%%%%

 
    Norm_form_sum = zero*zeros(4,4);
    % Beta is order 2 
    for i_beta = 0:2 
        j_beta = 2-i_beta ;
        A_beta = abs(normalForm_coeff(:,:,i_beta+1,j_beta+1));
        Norm_form_sum =Norm_form_sum +A_beta ;
    end
    Zc_vec = sum(Norm_form_sum) * K_Np2_u;
    Zc =max(Zc_vec)

    data_bndl_poly.Zc_vec =Zc_vec ;
 
    % Combine everything into the radii polynomial

    Y_0 = Ya_0+Yb_0+Yc_0
    try 
        Y_0.sup
        Y_0.inf
    end

    Z = Za+Zb+Zc

    data_bndl_poly.Y_0=Y_0;
    data_bndl_poly.Z=Z;

    if Z < 1
        r_min = Y_0 / (1-Z)
        if params.isIntval
            r_min = intval(r_min.sup)
        end
    else
        disp('Theorem failed! Z >= 1')
        r_min=nan;
    end

    data_bndl_poly.r_min = r_min;
    data_bndl_poly.K_N_u=K_N_u;
    data_bndl_poly.K_N_s=K_N_s;
end