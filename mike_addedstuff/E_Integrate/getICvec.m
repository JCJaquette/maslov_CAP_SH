function orthog = getICvec(pulseIC,basisVec1,basisVec2)
% project pulse IC into E^u_- and then find orthogonal vector in E^u_-

    A = [basisVec1,basisVec2];

    



    projectedIC = A*inv(A'*A)*A' * pulseIC;

    if dot(projectedIC,basisVec1) > 10^-3
        vec = basisVec1;
    else
        vec = basisVec2;
    end
  
    proj = dot(vec,projectedIC)/dot(projectedIC,projectedIC) * projectedIC;

    orthog = vec - proj;

    % dot(orthog,vec)
    % dot(orthog,proj)
    % dot(orthog,projectedIC)

end

