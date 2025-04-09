function orthog = getICvec(pulseIC,basisVec1,basisVec2)
% project pulse IC into E^u_- and then find orthogonal vector in E^u_-

    A = [basisVec1,basisVec2];

    projectedIC = A*inv(A'*A)*A' * pulseIC;

    % disp("Diff from phi ' and its projection")
    % pulseIC-projectedIC

    if dot(projectedIC,basisVec1) > 10^-3
        basisvec = basisVec1;
    else
        basisvec = basisVec2;
    end
  
    basis_proj = dot(basisvec,projectedIC)/dot(projectedIC,projectedIC) * projectedIC;

    orthog = basisvec - basis_proj;

    % dot(orthog,basis_proj)
    % dot(orthog,projectedIC)
    % dot(pulseIC,orthog)
    % x = A*inv(A'*A)*A' * orthog;
    % norm(x-orthog)


end

