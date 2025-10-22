% This function takes the derivative of vector field for the
% swift-hohenberg equation and plugs in the coefficients for the unstable
% manifold. 

function coeffs = DGPbundle(params, mflds)

    % Make for either stable / unstable ... eventually
    P = mflds.stable.coeffs;
   
    p1 = P(:,:,1);
    
    order = params.mfld.order;
    coeffs = zeros(4,4,order+1,order+1);
    
    if params.isIntval
        coeffs =intval(coeffs);
    end
    
    coeffs(:,:,1,1) = [0,1,0,0; 0,0,1,0; 0,0,0,1; -params.mu - 1, 0, -2, 0];

    for alpha = 1:order 
        for i = 0:alpha 
            j = alpha - i;  
            coeffs41 = 2*params.nu*p1(i+1,j+1) - 3*starVec(p1,p1,i,j);
            coeffs(:,:,i+1,j+1) = [0,0,0,0;0,0,0,0;0,0,0,0;coeffs41, 0, 0, 0];
        end
    end
    
end
