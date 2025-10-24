function K_alpha = K_op(e_val,m,n)
% K_alpha - but only multiplied against the coefficients, not the normal
% form. 
% Resonances are set to zero. 
    % e_val -- a row vector
    % e_val = [values.s, values.u];
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