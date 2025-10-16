function [ r_min ] = bundle_rad_poly(params,mflds,bndl)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
    All_Bundle_coeffs = bndl.coeffs;
    normalForm_coeff = bndl.normalForm;

    order = params.mfld.order; 
    G_hat = DFQbundle(params, mflds); % I think this is \hat G

    values = mflds.values;
    e_val = [values.s, values.u];

    order_2 = 3*(order+1)-1;

    params_extend = params;
    params_extend.mfld.order=order_2 ;

    if params.isIntval
        zero=intval(0);
    else
        zero=0;
    end

    mflds_extend =mflds;

    Q =  zeros(order_2 +1,order_2 +1,4)*zero;
    Q(1:order +1,1:order +1,:)=mflds.stable.coeffs;
    mflds_extend.stable.coeffs=Q;

    All_Bundle_coeffs_extend =  zeros(4,4,order_2 +1,order_2 +1)*zero;
    All_Bundle_coeffs_extend(:,:,1:order +1,1:order +1) =All_Bundle_coeffs;
 

    G_hat_N = DFQbundle(params_extend, mflds_extend); % I think this is \hat G
    
    p1 = mflds.stable.coeffs(:,:,1);


    % W_0 
    W_0 = All_Bundle_coeffs(:,:,1,1);
    W_0_inv = inv(W_0);

    W_0_norm = norm(W_0 ,1);
    W_0_inv_norm = norm(W_0_inv ,1);

    % P_infty
    P_infty = + mflds.stable.r_min;


    % Epsilon_infty
    % Lemma 6.7

    eps_infty = 2*params.nu* P_infty +3*(2*sum(abs(p1),'all')*P_infty +P_infty^2 );
    
    % K_N 
    % Eq (3.7)
    [rho, theta]= rTh_coord(params) ;


    K_N = 1 / ( order * sqrt(rho) * cos( theta/2) );
    K_Np2 = 1 / ( (order +2)* sqrt(rho) * cos( theta/2) );

%%% More bounds 
G_N_norm = sum(abs(G_hat_N),'all'); 
W_N_norm = sum(abs(All_Bundle_coeffs),'all'); 
A_norm = sum(abs(normalForm_coeff),'all'); 

    %%%  Y bound  %%% 

    %%%  Y_0^a  %%% 

    % little_sum = zero; 
    little_sum_2 = zero; 
    for alpha = order+1:3*order 
        for i = 0:alpha 
            j = alpha - i;  
            % summand = starMat(G_hat_N,All_Bundle_coeffs_extend,i,j);
            summand_2 = K_op(e_val,i,j).* (W_0\ starMat(G_hat_N,All_Bundle_coeffs_extend,i,j));
            % little_sum =little_sum +sum(abs(summand),'all'); 
            little_sum_2 =little_sum_2 +norm(abs(summand_2),1); 
        end
    end





    % This needs to compute the convolution and throw away lower order terms
    % Ya_0 = K_N *W_0_inv_norm  *little_sum
    Ya_0 = little_sum_2



    %%%  Y_0^b  %%% 

    Yb_0 = K_N*W_0_inv_norm * eps_infty* W_N_norm

    %%%  Y_0^c  %%% 
% This only needs things of order N+1 and N+2
    
    little_sum = zero; 
    for alpha = order-1:order 
        for i = 0:alpha 
            j = alpha - i;  
            little_sum =little_sum +sum(abs(All_Bundle_coeffs(:,:,i+1,j+1)),'all');
        end
    end
    
    Yc_0= K_N * A_norm *little_sum

G_hat_N_norm=zeros(order_2 +1,order_2 +1);
    for i = 1:order_2 +1
        for j= 1:order_2 +1
            G_hat_N_norm(i,j) = norm(G_hat_N(:,:,i,j),1);
        end
    end

    surf(G_hat_N_norm)

    %%%  Z bound  %%% 

    %%%  Z_0^a  %%% 

    % first_component = sum(abs(G_hat_N(:,:,1,1)),'all');

    little_sum = zero; 
    for alpha = 1:order_2 
        % K_N_alpha = 1 / (( order+alpha) * sqrt(rho) * cos( theta/2) ) ;
        for i = 0:alpha 
            j = alpha - i;  
            % summand = starMat(G_hat_N,All_Bundle_coeffs_extend,i,j);
            summand =  (W_0\ G_hat_N(:,:,i+1,j+1)*W_0 );
            % little_sum =little_sum +sum(abs(summand),'all'); 
            little_sum =little_sum +norm(abs(summand),1); 
            % little_sum_2 =little_sum_2 +norm(abs(summand),1)*K_N_alpha ; 
        end
    end

 

    % % % % Za =  (G_N_norm-first_component )*K_N*W_0_inv_norm*W_0_norm
    Za =  little_sum*K_N
    % Za =  little_sum_2 


    %%%  Z_0^b  %%% 

    Zb =  K_N*eps_infty*W_0_inv_norm*W_0_norm

    %%%  Z_0^c  %%% 
    Zc =  A_norm*K_Np2

    %  Rad Poly 

    Y_0 = Ya_0+Yb_0+Yc_0

    Z = Za+Zb+Zc

    if Z < 1
        r_min = Y_0 / (1-Z)
    else
        disp('Theorem failed! Z >= 1')
        r_min=nan;
    end

end