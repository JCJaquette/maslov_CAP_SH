function [rho, theta]= rTh_coord(params) 
    % Stores 

    % Eq (3.7) 
    rho = (1+abs(params.mu))^(1/2);    


    if params.isIntval
        theta = atan(-sqrt(params.mu))+intval('pi');
        % theta2 = atan2(sqrt(params.mu).mid,-1);
    else
        theta = atan2(sqrt(params.mu),-1);
        % theta2 = atan(-sqrt(params.mu))+intval('pi');
    end
 
    if theta <pi/2 || theta >pi
        disp('ERROR ERROR: Wrong branch of arctan')
    end
    
end