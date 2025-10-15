function [ r_min ] = bundle_rad_poly(params,mflds,bndl)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
    All_Bundle_coeffs = bndl.coeffs;
    normalForm_coeff = bndl.normalForm;

    order = params.mfld.order; 
    G_hat = DFQbundle(params, mflds); % I think this is \hat G

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

    W_0_norm = norm(W_0 );
    W_0_inv_norm = norm(W_0_inv );

    % P_infty
    P_infty = + mflds.stable.r_min;


    % Epsilon_infty
    % Lemma 6.7

    eps_infty = 2*params.nu* P_infty +3*(2*sum(abs(p1),'all')*P_infty +P_infty^2 );
    
    % K_N 
    % Eq (3.7)
    [rho, theta]= rTh_coord(params) ;


    K_N = 1 / ( order + sqrt(rho) * cos( theta/2) );

%%% More bounds 
G_N_norm = sum(abs(G_hat_N),'all'); 
W_N_norm = sum(abs(All_Bundle_coeffs),'all'); 
A_norm = sum(abs(normalForm_coeff),'all'); 

    %%%  Y bound  %%% 

    %%%  Y_0^a  %%% 

    little_sum = zero; 
    for alpha = order+1:3*order 
        for i = 0:alpha 
            j = alpha - i;  
            summand = starMat(G_hat_N,All_Bundle_coeffs_extend,i,j);
            little_sum =little_sum +sum(abs(summand),'all'); 
        end
    end

    % This needs to compute the convolution and throw away lower order terms
    Ya_0 = K_N *W_0_inv_norm  *little_sum

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


    %%%  Z bound  %%% 

    %%%  Z_0^a  %%% 
    first_component = sum(abs(G_hat_N(:,:,1,1)),'all');

    Za =  (G_N_norm-first_component )*K_N*W_0_inv_norm*W_0_norm
% Note: W_0_inv_norm*W_0_norm ~ 9, so multiply this out, and 
% subtract off the first coeffieinet 
% This'll get us nearly 2 orders of magnitude

    %%%  Z_0^b  %%% 

    Zb =  K_N*eps_infty*W_0_inv_norm*W_0_norm

    %%%  Z_0^c  %%% 
    Zc =  A_norm*K_N

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