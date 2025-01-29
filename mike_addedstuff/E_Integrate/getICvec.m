function orthog = getICvec(pulseIC,basisVec1,basisVec2)
% project pulse IC into E^u_- and then find orthogonal vector

    A = [basisVec1,basisVec2];
    projectedIC = A*inv(A'*A)*A' * pulseIC;

    if dot(projectedIC,basisVec1) > 10^-3
        vec = basisVec1;
    else
        vec = basisVec2;
    end
  
    proj = dot(projectedIC,vec)/dot(vec,vec) * vec;

    orthog = projectedIC - proj;

end

