% NEWTON_HALFLINE  Perform Newton's method on the pulse solution
% coefficients on domain [0,L]
%   S = S.Newton_halfline()
%   S = Newton)_halfline(S) 
% 
function full_uout = Newton_halfline(S) 
    % Option to display output and use Jacobian
    options=optimset('Display','iter','Jacobian','on','MaxIter',10000);     
    
    % Sometimes need to tweak the scaling of the normal form solution for
    % convergence to a pulse 

    u = S.normalForm.sol(:,1); 
    [uout,fval] = fsolve(@(u) fourierODE_halfline(S,u),u,options);  
    full_uout = [flip(uout); uout(2:end)];
    
end