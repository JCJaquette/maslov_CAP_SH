function [r theta]= Compute_Matrices(params,mflds) 
    % Stores 

    % Eq (3.7) 
    r = (1+abs(params.mu))^(1/2);    

    theta = atan2(sqrt(params.mu),-1);
    if theta <pi/2 || theta >pi
        disp('ERROR ERROR: Wrong branch of arctan')
    end
    
end